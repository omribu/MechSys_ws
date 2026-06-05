"""Parts 5 & 6: discrete-time (indirect digital) MRAC with ZOH; sampling-time
degradation and a PID phase-lead augmentation to recover performance."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sim_core import (plant_step, regressors, ref_step, A1, A2, A3, B3, G, C, H,
                      LAM, L1, L2, L3, M, B, K, ALPHA, TAU)

plt.rcParams.update({"figure.dpi": 120, "font.size": 11, "axes.grid": True,
                     "grid.alpha": 0.3})

DT_PLANT = 1e-3          # continuous plant integration step
T = 12.0
t = np.arange(0.0, T, DT_PLANT)
N = len(t)
r = np.ones(N)

def run_discrete(Ts, gamma=5.0, use_pid=False,
                 Kp=0.8, Ki=0.1, Kd=6.0, tau_d=0.05):
    """
    Indirect digital design: the continuous MRAC law is discretized. Output
    derivatives are obtained with discrete-time *filtered* differentiators (dirty
    derivatives, first-order low-pass pole at 1/tau_d) instead of raw differences,
    the standard robust choice for digital control of higher-relative-degree
    plants. Control is held by a ZOH over each period Ts while the continuous
    plant integrates at DT_PLANT. Only a_i, b_i are unknown; h = 1/c is known.
    """
    steps_per_sample = max(1, int(round(Ts/DT_PLANT)))
    state = np.array([0.0,0.0,0.0])
    ym_st = np.array([0.0,0.0,0.0])
    th_hat = np.zeros(5)
    u_held = 0.0
    integ = 0.0
    y_prev = 0.0; yd_f = 0.0; ydd_f = 0.0; yd_prev = 0.0
    y_hist=np.zeros(N); ym_hist=np.zeros(N); e_hist=np.zeros(N); u_hist=np.zeros(N)
    th_hist=np.zeros((N,5))
    a_f = tau_d/(tau_d+Ts)   # discrete low-pass coefficient

    for i in range(N):
        # ----- sample & recompute control at sampling instants -----
        if i % steps_per_sample == 0:
            x1,x2,x3 = state
            y = x1
            ym, ymd, ymdd = ym_st
            ymddd = L3*r[i] - L1*ymdd - L2*ymd - L3*ym
            # filtered derivatives (dirty differentiation)
            raw_yd = (y - y_prev)/Ts
            yd_f  = a_f*yd_f + (1-a_f)*raw_yd
            raw_ydd = (yd_f - yd_prev)/Ts
            ydd_f = a_f*ydd_f + (1-a_f)*raw_ydd
            yd, ydd = yd_f, ydd_f
            e   = y - ym
            ed  = yd - ymd
            edd = ydd - ymdd
            phi = regressors(y, yd, ydd)
            v   = ymddd - LAM[0]*edd - LAM[1]*ed - LAM[2]*e
            z   = edd + LAM[0]*ed + LAM[1]*e
            u_mrac = H*v + th_hat @ phi          # h = H known
            if use_pid:
                integ += e*Ts
                u_pid = -(Kp*e + Ki*integ + Kd*ed)   # error e=y-ym, drive to 0
                u_held = u_mrac + u_pid
            else:
                u_held = u_mrac
            # adaptation (discrete update) with sigma-normalization
            norm = 1.0 + phi@phi
            th_hat = th_hat - gamma*z*phi*Ts/norm
            y_prev = y; yd_prev = yd_f
        # ----- ZOH: hold u_held, integrate continuous plant -----
        u_hist[i] = u_held
        y_hist[i] = state[0]
        ym_hist[i] = ym_st[0]
        e_hist[i] = state[0] - ym_st[0]
        th_hist[i] = th_hat
        state = plant_step(state, u_held, DT_PLANT)
        ym_st, _ = ref_step(ym_st, r[i], DT_PLANT)
    return dict(y=y_hist, ym=ym_hist, e=e_hist, u=u_hist, th=th_hist)


# also need the continuous adaptive baseline for overlay (h known)
def run_continuous(gamma=5.0):
    state=np.array([0.0,0.0,0.0]); ym_st=np.array([0.0,0.0,0.0])
    th_hat=np.zeros(5)
    y=np.zeros(N); ym=np.zeros(N); e=np.zeros(N); u=np.zeros(N)
    for i in range(N):
        x1,x2,x3=state
        yv=x1; ydv=x2
        yddv=(1.0/M)*(x3-B*x2-K*x1-ALPHA*x1**3)
        m0,m1,m2=ym_st
        ymddd=L3*r[i]-L1*m2-L2*m1-L3*m0
        ev=yv-m0; edv=ydv-m1; eddv=yddv-m2
        z=eddv+LAM[0]*edv+LAM[1]*ev
        phi=regressors(yv,ydv,yddv)
        v=ymddd-LAM[0]*eddv-LAM[1]*edv-LAM[2]*ev
        uu=H*v+th_hat@phi
        u[i]=uu; y[i]=yv; ym[i]=m0; e[i]=ev
        th_hat=th_hat-gamma*z*phi*DT_PLANT
        state=plant_step(state,uu,DT_PLANT)
        ym_st,_=ref_step(ym_st,r[i],DT_PLANT)
    return dict(y=y,ym=ym,e=e,u=u)

if __name__ == "__main__":
    cont = run_continuous()
    Ts20 = 1/20.0
    dis = run_discrete(Ts20)

    # Fig P5: discrete control signal (ZOH staircase)
    plt.figure(figsize=(7,4))
    plt.plot(t, dis['u'], 'C0', label='u (discrete + ZOH, 20 Hz)')
    plt.xlabel('t [s]'); plt.ylabel('u(t)'); plt.title('Part 5: discrete control signal (ZOH)')
    plt.xlim(0,4); plt.legend(); plt.tight_layout(); plt.savefig('figs/p5_u_zoh.png'); plt.close()

    # Fig P5: discrete vs continuous tracking
    plt.figure(figsize=(7,4))
    plt.plot(t, dis['ym'], 'C0', label='reference model')
    plt.plot(t, cont['y'], 'C2', label='continuous plant')
    plt.plot(t, dis['y'],  'C1', label='discrete plant (20 Hz)')
    plt.plot(t, cont['e'], 'C5', lw=0.8, label='continuous error')
    plt.plot(t, dis['e'],  'C3', lw=0.8, label='discrete error')
    plt.xlabel('t [s]'); plt.ylabel('amplitude'); plt.title('Part 5: discrete vs continuous controller')
    plt.legend(fontsize=8); plt.tight_layout(); plt.savefig('figs/p5_compare.png'); plt.close()

    # Fig P5: discrete parameter estimates
    plt.figure(figsize=(7,4))
    labels=['$\\hat a_1$','$\\hat a_2$','$\\hat a_3$','$\\hat b_3$','$\\hat g$']
    for j in range(5): plt.plot(t, dis['th'][:,j], label=labels[j])
    plt.xlabel('t [s]'); plt.ylabel('estimate'); plt.title('Part 5: parameter estimates (discrete)')
    plt.legend(fontsize=8,ncol=2); plt.tight_layout(); plt.savefig('figs/p5_params.png'); plt.close()

    # Fig P6: sampling-time sweep (two panels: tracking zoom + RMSE vs Ts)
    Tss=[0.02,0.035,0.05,0.06,0.065]
    def _rmse(e):
        e=np.asarray(e); fin=np.isfinite(e)
        return float(np.sqrt(np.mean(e[fin]**2))) if fin.any() else float('nan')
    fig,(axA,axB)=plt.subplots(1,2,figsize=(10,4))
    axA.plot(t, cont['ym'],'k',lw=2.0,label='reference model')
    rmses=[]
    for Ts in Tss:
        d=run_discrete(Ts)
        axA.plot(t, d['y'], lw=1.0, label=f'{Ts*1000:.0f} ms ({1/Ts:.0f} Hz)')
        rmses.append(_rmse(d['e']))
    axA.set_xlim(0,3.5); axA.set_ylim(0,1.15)
    axA.set_xlabel('t [s]'); axA.set_ylabel('amplitude')
    axA.set_title('Tracking (transient zoom)'); axA.legend(fontsize=8)
    axB.bar([f'{Ts*1000:.0f}' for Ts in Tss], rmses, color='C3', alpha=0.8)
    axB.set_xlabel('sampling time Ts [ms]'); axB.set_ylabel('tracking RMSE')
    axB.set_title('Performance degradation vs Ts')
    for i,v in enumerate(rmses):
        axB.text(i, v, f'{v:.3f}', ha='center', va='bottom', fontsize=8)
    plt.tight_layout(); plt.savefig('figs/p6_sweep.png'); plt.close()

    # Fig P6: PID-augmented vs plain discrete at a coarse rate (50 ms)
    Ts_coarse = 0.05
    d_plain = run_discrete(Ts_coarse, use_pid=False)
    d_pid   = run_discrete(Ts_coarse, use_pid=True)
    plt.figure(figsize=(7,4))
    plt.plot(t, d_plain['ym'],'k',lw=1.4,label='reference model')
    plt.plot(t, d_plain['y'],'C1',label='discrete MRAC')
    plt.plot(t, d_pid['y'],'C2',label='discrete MRAC + PID')
    plt.plot(t, d_plain['e'],'C3',lw=0.8,label='error (MRAC)')
    plt.plot(t, d_pid['e'],'C0',lw=0.8,label='error (MRAC+PID)')
    plt.xlabel('t [s]'); plt.ylabel('amplitude'); plt.title('Part 6: PID augmentation recovers performance (Ts=50 ms)')
    plt.xlim(0,6); plt.legend(fontsize=8); plt.tight_layout(); plt.savefig('figs/p6_pid.png'); plt.close()

    # Fig P6: control signals plain vs PID
    plt.figure(figsize=(7,4))
    plt.plot(t, d_plain['u'],'C1',label='u (MRAC)')
    plt.plot(t, d_pid['u'],'C2',label='u (MRAC+PID)')
    plt.xlabel('t [s]'); plt.ylabel('u(t)'); plt.title('Part 6: control signal with/without PID')
    plt.xlim(0,4); plt.legend(); plt.tight_layout(); plt.savefig('figs/p6_pid_u.png'); plt.close()

    # Fig P6: parameter estimates plain vs PID
    plt.figure(figsize=(7,4))
    labels=['$\\hat a_1$','$\\hat a_2$','$\\hat a_3$','$\\hat b_3$','$\\hat g$']
    cols=['C0','C1','C2','C3','C4']
    for j in range(5):
        plt.plot(t, d_plain['th'][:,j], cols[j], lw=1.0, label=labels[j]+' (MRAC)')
        plt.plot(t, d_pid['th'][:,j], cols[j], lw=1.0, ls='--', label=labels[j]+' (+PID)')
    plt.xlabel('t [s]'); plt.ylabel('estimate'); plt.title('Part 6: coefficient estimates with/without PID')
    plt.legend(fontsize=7, ncol=2); plt.tight_layout(); plt.savefig('figs/p6_pid_params.png'); plt.close()

    import numpy as _np
    def rmse(e):
        e=_np.asarray(e); fin=_np.isfinite(e)
        return float(_np.sqrt(_np.mean(e[fin]**2))) if fin.any() else float('nan')
    print("continuous RMSE:", round(rmse(cont['e']),5))
    print("discrete 20Hz RMSE:", round(rmse(dis['e']),5))
    print("discrete 50ms plain RMSE:", round(rmse(d_plain['e']),5))
    print("discrete 50ms +PID RMSE:", round(rmse(d_pid['e']),5))
    for Ts in Tss:
        d=run_discrete(Ts); print(f"  Ts={Ts*1000:.0f}ms RMSE={rmse(d['e']):.5f}")
    print("Figures written.")
