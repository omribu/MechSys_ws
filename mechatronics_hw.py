"""
mechatronics_hw.py
==================
Nonlinear Adaptive Control (MRAC) of a third-order mass-spring-damper system
with a nonlinear (cubic) spring and a first-order actuator.

Notation matches the report exactly:
    plant (normalized):  y''' + a1 y'' + a2 y' + a3 y + b1 (y^2 y') + b2 y^3 = c u
    parameter vector  theta = [a1, a2, a3, b1, b2]
    regressor vector  phi   = [y'', y', y, y^2 y', y^3]
    filtered error    s = e'' + alpha2 e' + alpha1 e,    e = y - y_m
    controller gains  k2 = alpha2 + lam,  k1 = alpha1 + lam alpha2,  k0 = lam alpha1
    adaptation        theta_hat' = -Gamma phi s

Run:  python3 mechatronics_hw.py
Generates every figure used in the report into ./figures/
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"figure.dpi": 120, "font.size": 11,
                     "axes.grid": True, "grid.alpha": 0.3})

# ======================================================================
# 1. PHYSICAL PARAMETERS AND NORMALIZED MODEL COEFFICIENTS
# ======================================================================
m, b, k, alpha, tau = 1.0, 0.5, 2.0, 1.0, 0.5

a1 = (tau * b + m) / (tau * m)     # = 2.5   coeff of y''
a2 = (tau * k + b) / (tau * m)     # = 3.0   coeff of y'
a3 = k / (tau * m)                 # = 4.0   coeff of y
b1 = 3.0 * alpha / m               # = 3.0   coeff of y^2 y'  (mixed nonlinearity)
b2 = alpha / (tau * m)             # = 2.0   coeff of y^3     (cubic nonlinearity)
c  = 1.0 / (tau * m)               # = 2.0   input gain

THETA = np.array([a1, a2, a3, b1, b2])   # true parameter vector


def phi_vec(y, yd, ydd):
    """Regressor phi = [y'', y', y, y^2 y', y^3]."""
    return np.array([ydd, yd, y, y * y * yd, y ** 3])


# ======================================================================
# 2. REFERENCE MODEL  (3rd order, Hurwitz, unity DC gain)
#    y_m''' + lam2 y_m'' + lam1 y_m' + lam0 y_m = lam0 r
#    chosen as a triple real pole at p:  (s + p)^3
# ======================================================================
p_ref = 2.5
lam2 = 3.0 * p_ref            # 7.5
lam1 = 3.0 * p_ref ** 2       # 18.75
lam0 = p_ref ** 3             # 15.625

# ======================================================================
# 3. FILTERED-ERROR AND CONTROLLER GAINS
#    s = e'' + alpha2 e' + alpha1 e   (Hurwitz: p^2 + alpha2 p + alpha1)
#    k2 = alpha2 + lam, k1 = alpha1 + lam alpha2, k0 = lam alpha1
# ======================================================================
alpha1 = 6.0
alpha2 = 5.0
lam    = 5.0
k2 = alpha2 + lam
k1 = alpha1 + lam * alpha2
k0 = lam * alpha1


# ======================================================================
# CONTINUOUS-TIME CLOSED LOOP (Parts 3 & 4)
# ======================================================================
def simulate_continuous(T=12.0, dt=1e-3, adaptive=False, gamma=5.0, r_amp=1.0):
    """
    Integrates plant + reference model + controller with explicit Euler.
    Plant states:     x1=y, x2=y', x3=y''  (normalized companion form)
    Reference states: z1=y_m, z2=y_m', z3=y_m''
    If adaptive=False -> nominal control with known THETA.
    If adaptive=True  -> theta_hat adapted on-line (theta_hat(0)=0).
    """
    n = int(T / dt)
    t = np.arange(n) * dt
    r = r_amp * np.ones(n)

    x1 = x2 = x3 = 0.0          # plant
    z1 = z2 = z3 = 0.0          # reference model
    th = THETA.copy() if not adaptive else np.zeros(5)

    Y = np.zeros(n); YM = np.zeros(n); E = np.zeros(n); U = np.zeros(n)
    TH = np.zeros((n, 5))

    for i in range(n):
        # reference-model derivative
        z3dot = lam0 * r[i] - lam2 * z3 - lam1 * z2 - lam0 * z1
        # tracking error and its derivatives
        e   = x1 - z1
        ed  = x2 - z2
        edd = x3 - z3
        s   = edd + alpha2 * ed + alpha1 * e
        phi = phi_vec(x1, x2, x3)
        # control law:  u = (1/c)[ y_m''' - k2 e'' - k1 e' - k0 e + theta_hat^T phi ]
        u = (1.0 / c) * (z3dot - k2 * edd - k1 * ed - k0 * e + th @ phi)

        Y[i], YM[i], E[i], U[i], TH[i] = x1, z1, e, u, th

        # plant derivative (true parameters):  y''' = c u - THETA^T phi
        x3dot = c * u - THETA @ phi
        # integrate plant
        x1 += x2 * dt; x2 += x3 * dt; x3 += x3dot * dt
        # integrate reference model
        z1 += z2 * dt; z2 += z3 * dt; z3 += z3dot * dt
        # adaptation:  theta_hat' = -Gamma phi s
        if adaptive:
            th = th - gamma * phi * s * dt

    return dict(t=t, y=Y, ym=YM, e=E, u=U, th=TH)


def simulate_openloop(T=12.0, dt=1e-3, u_const=1.0):
    n = int(T / dt); t = np.arange(n) * dt
    x1 = x2 = x3 = 0.0; Y = np.zeros(n)
    for i in range(n):
        Y[i] = x1
        x3dot = c * u_const - THETA @ phi_vec(x1, x2, x3)
        x1 += x2 * dt; x2 += x3 * dt; x3 += x3dot * dt
    return t, Y


# ======================================================================
# DISCRETE-TIME CLOSED LOOP (Parts 5 & 6) -- indirect digital design
#   Continuous law discretized; derivatives via filtered (dirty) differences;
#   control held by ZOH; continuous plant integrated at dt_plant.
# ======================================================================
def simulate_discrete(Ts, T=12.0, dt_plant=1e-3, adaptive=True, gamma=5.0,
                      tau_d=0.03, use_pid=False, Kp=0.8, Ki=0.1, Kd=3.0,
                      r_amp=1.0):
    n = int(T / dt_plant); t = np.arange(n) * dt_plant
    r = r_amp * np.ones(n)
    steps = max(1, int(round(Ts / dt_plant)))

    x1 = x2 = x3 = 0.0
    z1 = z2 = z3 = 0.0
    th = THETA.copy() if not adaptive else np.zeros(5)
    u_hold = 0.0
    integ = 0.0

    # filtered-derivative memory (controller only sees sampled y)
    y_prev = 0.0; yd_f = 0.0; ydd_f = 0.0; yd_prev = 0.0
    af = tau_d / (tau_d + Ts)

    Y = np.zeros(n); YM = np.zeros(n); E = np.zeros(n); U = np.zeros(n)
    TH = np.zeros((n, 5))

    for i in range(n):
        if i % steps == 0:                      # ---- sampling instant ----
            y = x1
            z3dot = lam0 * r[i] - lam2 * z3 - lam1 * z2 - lam0 * z1
            # filtered (dirty) derivatives of the measured output
            raw_yd = (y - y_prev) / Ts
            yd_f = af * yd_f + (1 - af) * raw_yd
            raw_ydd = (yd_f - yd_prev) / Ts
            ydd_f = af * ydd_f + (1 - af) * raw_ydd
            yd, ydd = yd_f, ydd_f
            e   = y - z1
            ed  = yd - z2
            edd = ydd - z3
            s   = edd + alpha2 * ed + alpha1 * e
            phi = phi_vec(y, yd, ydd)
            u_mrac = (1.0 / c) * (z3dot - k2 * edd - k1 * ed - k0 * e + th @ phi)
            if use_pid:
                integ += e * Ts
                u_pid = -(Kp * e + Ki * integ + Kd * ed)
                u_hold = u_mrac + u_pid
            else:
                u_hold = u_mrac
            if adaptive:                        # discrete adaptation (normalized)
                norm = 1.0 + phi @ phi
                th = th - gamma * phi * s * Ts / norm
            y_prev = y; yd_prev = yd_f

        # ---- ZOH: hold u_hold, integrate continuous plant & reference ----
        Y[i], YM[i], E[i], U[i], TH[i] = x1, z1, x1 - z1, u_hold, th
        x3dot = c * u_hold - THETA @ phi_vec(x1, x2, x3)
        x1 += x2 * dt_plant; x2 += x3 * dt_plant; x3 += x3dot * dt_plant
        z3dot2 = lam0 * r[i] - lam2 * z3 - lam1 * z2 - lam0 * z1
        z1 += z2 * dt_plant; z2 += z3 * dt_plant; z3 += z3dot2 * dt_plant

    return dict(t=t, y=Y, ym=YM, e=E, u=U, th=TH)


def rmse(e):
    e = np.asarray(e); fin = np.isfinite(e)
    return float(np.sqrt(np.mean(e[fin] ** 2))) if fin.any() else float("nan")


# ======================================================================
# MAIN: run everything and produce the figures
# ======================================================================
if __name__ == "__main__":
    LBL = ['$\\hat a_1$', '$\\hat a_2$', '$\\hat a_3$', '$\\hat b_1$', '$\\hat b_2$']

    # ---------- Part 3: nominal (known parameters) ----------
    nom = simulate_continuous(adaptive=False)
    t_ol, y_ol = simulate_openloop()

    plt.figure(figsize=(7, 4))
    plt.plot(t_ol, np.ones_like(t_ol), 'C0', label='reference (step)')
    plt.plot(t_ol, y_ol, 'C1', label='plant (open loop)')
    plt.xlabel('t [s]'); plt.ylabel('amplitude')
    plt.title('Open-loop step response (no controller)')
    plt.legend(); plt.tight_layout()
    plt.savefig('figures/openloop.png'); plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(nom['t'], nom['ym'], 'C0', label='reference model $y_m$')
    plt.plot(nom['t'], nom['y'], 'C2', label='plant $y$ (known params)')
    plt.plot(nom['t'], nom['e'], 'C3', label='tracking error $e$')
    plt.xlabel('t [s]'); plt.ylabel('amplitude')
    plt.title('Nominal control with known parameters')
    plt.legend(); plt.tight_layout()
    plt.savefig('figures/nominal_tracking.png'); plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(nom['t'], nom['u'], 'C4')
    plt.xlabel('t [s]'); plt.ylabel('$u(t)$')
    plt.title('Nominal control signal')
    plt.tight_layout(); plt.savefig('figures/nominal_u.png'); plt.close()

    # ---------- Part 4: adaptive (unknown parameters) ----------
    ada = simulate_continuous(adaptive=True, gamma=5.0)

    plt.figure(figsize=(7, 4))
    plt.plot(ada['t'], ada['ym'], 'C0', label='reference model $y_m$')
    plt.plot(nom['t'], nom['y'], 'C2', label='plant (known)')
    plt.plot(ada['t'], ada['y'], 'C1', label='plant (adaptive)')
    plt.plot(ada['t'], ada['e'], 'C3', lw=0.9, label='error (adaptive)')
    plt.xlabel('t [s]'); plt.ylabel('amplitude')
    plt.title('Adaptive vs. known-parameter tracking')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig('figures/adaptive_tracking.png'); plt.close()

    plt.figure(figsize=(7, 4))
    for j in range(5):
        plt.plot(ada['t'], ada['th'][:, j], label=LBL[j])
    plt.xlabel('t [s]'); plt.ylabel('estimate')
    plt.title('Adaptive parameter estimates')
    plt.legend(fontsize=8, ncol=2); plt.tight_layout()
    plt.savefig('figures/adaptive_params.png'); plt.close()

    # ---------- Part 5: discrete controller ----------
    cont = simulate_continuous(adaptive=True, gamma=5.0)
    dis = simulate_discrete(Ts=1 / 20.0, adaptive=True, gamma=5.0)

    plt.figure(figsize=(7, 4))
    plt.plot(dis['t'], dis['u'], 'C0', label='u (discrete + ZOH, 20 Hz)')
    plt.xlabel('t [s]'); plt.ylabel('$u(t)$'); plt.xlim(0, 4)
    plt.title('Discrete control signal (ZOH staircase)')
    plt.legend(); plt.tight_layout(); plt.savefig('figures/discrete_u.png'); plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(dis['t'], dis['ym'], 'C0', label='reference model')
    plt.plot(cont['t'], cont['y'], 'C2', label='continuous plant')
    plt.plot(dis['t'], dis['y'], 'C1', label='discrete plant (20 Hz)')
    plt.plot(dis['t'], dis['e'], 'C3', lw=0.8, label='discrete error')
    plt.xlabel('t [s]'); plt.ylabel('amplitude')
    plt.title('Discrete vs. continuous controller')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig('figures/discrete_tracking.png'); plt.close()

    plt.figure(figsize=(7, 4))
    for j in range(5):
        plt.plot(dis['t'], dis['th'][:, j], label=LBL[j])
    plt.xlabel('t [s]'); plt.ylabel('estimate')
    plt.title('Parameter estimates (discrete controller)')
    plt.legend(fontsize=8, ncol=2); plt.tight_layout()
    plt.savefig('figures/discrete_params.png'); plt.close()

    # ---------- Part 6: sampling-time sweep + PID ----------
    Tss = [0.02, 0.035, 0.05, 0.06, 0.065, 0.07]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10, 4))
    axA.plot(cont['t'], cont['ym'], 'k', lw=2.0, label='reference model')
    rmses = []
    for Ts in Tss:
        d = simulate_discrete(Ts=Ts, adaptive=True, gamma=5.0)
        rr = rmse(d['e'])
        rmses.append(rr)
        if np.isfinite(d['y']).all() and rr < 1.0:
            axA.plot(d['t'], d['y'], lw=1.0, label=f'{Ts*1000:.0f} ms')
        else:
            axA.plot([], [], lw=1.0, label=f'{Ts*1000:.0f} ms (unstable)')
    axA.set_xlim(0, 3.5); axA.set_ylim(0, 1.15)
    axA.set_xlabel('t [s]'); axA.set_ylabel('amplitude')
    axA.set_title('Tracking (transient zoom)'); axA.legend(fontsize=8)
    bar_vals = [rr if (np.isfinite(rr) and rr < 1.0) else 0.0 for rr in rmses]
    axB.bar([f'{Ts*1000:.0f}' for Ts in Tss], bar_vals, color='C3', alpha=0.8)
    axB.set_xlabel('sampling time $T_s$ [ms]'); axB.set_ylabel('tracking RMSE')
    axB.set_title('Performance degradation vs $T_s$')
    for ii, (rr, bv) in enumerate(zip(rmses, bar_vals)):
        ok = np.isfinite(rr) and rr < 1.0
        axB.text(ii, bv if ok else max(bar_vals) * 0.5,
                 f'{rr:.3f}' if ok else 'unstable',
                 ha='center', va='bottom', fontsize=8,
                 rotation=0 if ok else 90, color='black' if ok else 'C3')
    plt.tight_layout(); plt.savefig('figures/sampling_sweep.png'); plt.close()

    Ts_c = 0.05
    d_plain = simulate_discrete(Ts=Ts_c, adaptive=True, use_pid=False)
    d_pid   = simulate_discrete(Ts=Ts_c, adaptive=True, use_pid=True)

    plt.figure(figsize=(7, 4))
    plt.plot(d_plain['t'], d_plain['ym'], 'k', lw=1.4, label='reference model')
    plt.plot(d_plain['t'], d_plain['y'], 'C1', label='discrete MRAC')
    plt.plot(d_pid['t'], d_pid['y'], 'C2', label='discrete MRAC + PID')
    plt.plot(d_plain['t'], d_plain['e'], 'C3', lw=0.8, label='error (MRAC)')
    plt.plot(d_pid['t'], d_pid['e'], 'C0', lw=0.8, label='error (MRAC+PID)')
    plt.xlim(0, 6); plt.xlabel('t [s]'); plt.ylabel('amplitude')
    plt.title('PID augmentation recovers performance ($T_s=50$ ms)')
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig('figures/pid_tracking.png'); plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(d_plain['t'], d_plain['u'], 'C1', label='u (MRAC)')
    plt.plot(d_pid['t'], d_pid['u'], 'C2', label='u (MRAC+PID)')
    plt.xlim(0, 4); plt.xlabel('t [s]'); plt.ylabel('$u(t)$')
    plt.title('Control signal with / without PID')
    plt.legend(); plt.tight_layout(); plt.savefig('figures/pid_u.png'); plt.close()

    # ---------- console summary ----------
    print("Physical:  m=%.2f b=%.2f k=%.2f alpha=%.2f tau=%.2f" % (m, b, k, alpha, tau))
    print("Normalized theta = [a1,a2,a3,b1,b2] =", THETA, " c =", c)
    print("Reference poles: triple at p =", p_ref,
          " (lam2,lam1,lam0)=(%.3f,%.3f,%.3f)" % (lam2, lam1, lam0))
    print("Nominal   final error  = %.3e" % nom['e'][-1])
    print("Adaptive  final error  = %.3e" % ada['e'][-1])
    print("Continuous RMSE        = %.5f" % rmse(cont['e']))
    print("Discrete 20Hz RMSE     = %.5f" % rmse(dis['e']))
    print("Discrete 50ms RMSE     = %.5f" % rmse(d_plain['e']))
    print("Discrete 50ms+PID RMSE = %.5f" % rmse(d_pid['e']))
    for Ts, rr in zip(Tss, rmses):
        print("   Ts=%3.0f ms  RMSE=%.5f" % (Ts * 1000, rr))
    print("All figures written to ./figures/")
