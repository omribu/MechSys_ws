"""Parts 3 & 4: continuous-time MRAC, known parameters vs adaptive (unknown)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sim_core import (plant_step, regressors, ref_step, A1, A2, A3, B3, G, C, H,
                      LAM, L1, L2, L3, M, B, K, ALPHA, TAU)

plt.rcParams.update({"figure.dpi": 120, "font.size": 11, "axes.grid": True,
                     "grid.alpha": 0.3})

dt = 1e-3
T  = 12.0
t  = np.arange(0.0, T, dt)
N  = len(t)
r  = np.ones(N)           # step reference into the reference model

# true normalized parameter vector theta = [a1,a2,a3,b3,g] and h = 1/c
THETA = np.array([A1, A2, A3, B3, G])

# numerical derivative helper using stored history (backward differences)
def deriv(arr, i, dt):
    if i < 1: return 0.0
    return (arr[i] - arr[i-1]) / dt

# ============================================================
# PART 3: nominal (known parameters) controller
#   u = h*( ym''' - lam1 e'' - lam2 e' - lam3 e ) + theta . phi
#   where e = y - ym (tracking error). This cancels plant dynamics exactly,
#   leaving the error governed by Hurwitz polynomial -> e -> 0.
# ============================================================
def run_known():
    state = np.array([0.0, 0.0, 0.0])      # x1,x2,x3
    ym_st = np.array([0.0, 0.0, 0.0])
    y_hist  = np.zeros(N); yd_hist = np.zeros(N); ydd_hist = np.zeros(N)
    ym_hist = np.zeros(N); e_hist  = np.zeros(N); u_hist   = np.zeros(N)
    e = ed = edd = 0.0
    for i in range(N):
        x1, x2, x3 = state
        y, yd = x1, x2
        ydd = (1.0/M)*(x3 - B*x2 - K*x1 - ALPHA*x1**3)
        y_hist[i], yd_hist[i], ydd_hist[i] = y, yd, ydd
        ym, ymd, ymdd = ym_st
        ym_hist[i] = ym
        ymddd = L3*r[i] - L1*ymdd - L2*ymd - L3*ym
        e   = y - ym
        ed  = yd - ymd
        edd = ydd - ymdd
        e_hist[i] = e
        phi = regressors(y, yd, ydd)
        # desired y''' so that e''' + lam1 e'' + lam2 e' + lam3 e = 0
        ydddes = ymddd - LAM[0]*edd - LAM[1]*ed - LAM[2]*e
        u = H*ydddes + THETA @ phi
        u_hist[i] = u
        state = plant_step(state, u, dt)
        ym_st, _ = ref_step(ym_st, r[i], dt)
    return dict(y=y_hist, ym=ym_hist, e=e_hist, u=u_hist)

# ============================================================
# PART 4: adaptive controller (unknown a_i, b_i, h)
#   u = h_hat*(ym''' - lam.e) + theta_hat . phi
#   z = e'' + lam1 e' + lam2 e   (filtered tracking error, rel. deg handled)
#   adaptation:  theta_hat' = -gamma * z * phi ;  h_hat' = -gamma_h * z * v
# ============================================================
def run_adaptive(gamma=5.0):
    state = np.array([0.0, 0.0, 0.0])
    ym_st = np.array([0.0, 0.0, 0.0])
    th_hat = np.zeros(5)      # estimates of [a1,a2,a3,b3,g]; h known (= H)
    y_hist=np.zeros(N); ym_hist=np.zeros(N); e_hist=np.zeros(N); u_hist=np.zeros(N)
    th_hist=np.zeros((N,5))
    for i in range(N):
        x1, x2, x3 = state
        y, yd = x1, x2
        ydd = (1.0/M)*(x3 - B*x2 - K*x1 - ALPHA*x1**3)
        ym, ymd, ymdd = ym_st
        ym_hist[i] = ym
        ymddd = L3*r[i] - L1*ymdd - L2*ymd - L3*ym
        e   = y - ym
        ed  = yd - ymd
        edd = ydd - ymdd
        e_hist[i] = e; y_hist[i] = y
        # filtered error z = edd + lam1 ed + lam2 e  (Hurwitz combination, SPR)
        z = edd + LAM[0]*ed + LAM[1]*e
        phi = regressors(y, yd, ydd)
        v   = ymddd - LAM[0]*edd - LAM[1]*ed - LAM[2]*e
        u = H*v + th_hat @ phi          # h = H known
        u_hist[i] = u
        th_hist[i] = th_hat
        # adaptation law (Lyapunov/gradient based) for a_i, b_i only
        th_hat = th_hat - gamma*z*phi*dt
        state = plant_step(state, u, dt)
        ym_st, _ = ref_step(ym_st, r[i], dt)
    return dict(y=y_hist, ym=ym_hist, e=e_hist, u=u_hist, th=th_hist)

# -------------------- open loop response (no controller) --------------------
def run_openloop():
    state = np.array([0.0,0.0,0.0]); y=np.zeros(N)
    for i in range(N):
        y[i]=state[0]
        state=plant_step(state,1.0,dt)   # constant unit command
    return y

if __name__ == "__main__":
    ol = run_openloop()
    kn = run_known()
    ad = run_adaptive()

    # Fig: open loop
    plt.figure(figsize=(7,4))
    plt.plot(t, np.ones(N), 'C0', label='reference (step)')
    plt.plot(t, ol, 'C1', label='plant (open loop)')
    plt.xlabel('t [s]'); plt.ylabel('amplitude'); plt.title('Open-loop step response (no controller)')
    plt.legend(); plt.tight_layout(); plt.savefig('figs/p3_openloop.png'); plt.close()

    # Fig: known-parameter tracking
    plt.figure(figsize=(7,4))
    plt.plot(t, kn['ym'], 'C0', label='reference model $y_m$')
    plt.plot(t, kn['y'],  'C2', label='plant $y$ (known params)')
    plt.plot(t, kn['e'],  'C3', label='error $\\tilde y$')
    plt.xlabel('t [s]'); plt.ylabel('amplitude'); plt.title('Part 3: tracking with known parameters')
    plt.legend(); plt.tight_layout(); plt.savefig('figs/p3_known.png'); plt.close()

    # Fig: known control signal
    plt.figure(figsize=(7,4))
    plt.plot(t, kn['u'], 'C4'); plt.xlabel('t [s]'); plt.ylabel('u(t)')
    plt.title('Part 3: nominal control signal'); plt.tight_layout()
    plt.savefig('figs/p3_u.png'); plt.close()

    # Fig: adaptive vs known
    plt.figure(figsize=(7,4))
    plt.plot(t, ad['ym'], 'C0', label='reference model $y_m$')
    plt.plot(t, kn['y'],  'C2', label='plant (known)')
    plt.plot(t, ad['y'],  'C1', label='plant (adaptive)')
    plt.plot(t, ad['e'],  'C3', label='error (adaptive)')
    plt.plot(t, kn['e'],  'C5', lw=0.8, label='error (known)')
    plt.xlabel('t [s]'); plt.ylabel('amplitude'); plt.title('Part 4: adaptive vs known parameters')
    plt.legend(fontsize=8); plt.tight_layout(); plt.savefig('figs/p4_compare.png'); plt.close()

    # Fig: parameter estimates
    plt.figure(figsize=(7,4))
    labels = ['$\\hat a_1$','$\\hat a_2$','$\\hat a_3$','$\\hat b_3$','$\\hat g$']
    for j in range(5):
        plt.plot(t, ad['th'][:,j], label=labels[j])
    plt.xlabel('t [s]'); plt.ylabel('estimate'); plt.title('Part 4: adaptive parameter estimates')
    plt.legend(fontsize=8, ncol=2); plt.tight_layout(); plt.savefig('figs/p4_params.png'); plt.close()

    print("known: final error = %.4e" % kn['e'][-1])
    print("adaptive: final error = %.4e" % ad['e'][-1])
    print("adaptive estimates final:", np.round(ad['th'][-1],3))
    print("true theta:", np.round([A1,A2,A3,B3,G],3), "h=", H)
    print("Figures written.")
