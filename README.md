# Mechatronic Systems — Adaptive Control Homework (Python code)

Simulation code for the homework: adaptive control (MRAC) of a nonlinear
**third-order** mass–spring–damper with a nonlinear (cubic) spring and a
first-order actuator.

## The plant

Mechanical subsystem (second order):

    m*x'' + b*x' + k*x + alpha*x^3 = F(t)

Actuator (first order, time constant tau):

    tau*F' + F = u(t)

Eliminating the internal force F(t) gives the combined **third-order** plant:

    tau*m*x''' + (tau*b+m)*x'' + (tau*k+b)*x' + k*x
               + 3*tau*alpha*x^2*x' + alpha*x^3 = u(t)

Physical parameters used: m=1, b=0.5, k=2, alpha=1, tau=0.5
→ normalized: a1=1.25, a2=1.5, a3=2.0, b3=1.0 (cubic), g=1.5 (mixed term),
  c=2, h=1/c=0.5.

The controller treats h as known and adapts a1, a2, a3, b3, g
(per the assignment: "only a_i, b_i are unknown").

## Files

| file              | what it does                                                        |
|-------------------|---------------------------------------------------------------------|
| `sim_core.py`     | Plant model, regressors, reference model, shared constants.         |
| `sim_p34.py`      | **Part 3** (known params) and **Part 4** (adaptive) — continuous.   |
| `sim_p56.py`      | **Part 5** (discrete/indirect digital) and **Part 6** (Ts sweep + PID). |
| `make_diagrams.py`| Physical schematic, MRAC block diagram, reference-model step.       |
| `render_eqs.py`   | (Optional) renders the report equations to PNG — not needed to run sims. |
| `run_all.py`      | Runs everything and writes all figures to `./figs/`.                |

## How to run

Requirements: Python 3, `numpy`, `scipy`, `matplotlib` (and `sympy` only if you
want to re-verify the derivation; not required for the simulations).

    pip install numpy scipy matplotlib
    python3 run_all.py

All figures are written to `./figs/`. You can also run any part on its own:

    python3 sim_p34.py     # Parts 3 & 4
    python3 sim_p56.py     # Parts 5 & 6

## Key design choices (reflected in the code)

- **Reference model:** third order, triple real pole at p = 2.5 rad/s
  (lambda1 = 7.5, lambda2 = 18.75, lambda3 = 15.625). Fast enough for good
  tracking, slow enough to be compatible with the digital sampling rates.
- **Adaptation:** gradient/Lyapunov law  theta_hat' = -gamma * z * phi,
  with gamma = 5 and filtered error  z = e'' + lambda1 e' + lambda2 e.
- **Digital implementation (Part 5):** continuous law discretized (indirect
  design). Output derivatives use **filtered ("dirty") differentiators**
  (first-order low-pass, pole 1/tau_d) instead of raw differences — essential
  for stability of this relative-degree-3 loop. Control held by a ZOH.
- **Part 6 remedy:** a parallel PID (Kp=0.8, Ki=0.1, Kd=6.0) supplies phase
  lead that compensates the ZOH lag and restores performance at coarse Ts.

## What the numbers show

- Part 3 (known): tracking error ~1e-13 (exact cancellation).
- Part 4 (adaptive): error → ~1e-2 after a short learning transient.
  Estimates converge to tracking-optimal values, **not** the true parameters
  (a step input is not persistently exciting — expected MRAC behaviour).
- Part 5 (discrete, 20 Hz): RMSE ~0.037 vs ~0.021 continuous.
- Part 6: RMSE rises 0.031 → 0.042 from 20 ms to 60 ms, then jumps to 0.34 at
  65 ms (stability cliff). PID augmentation cuts the 50 ms RMSE 0.037 → 0.025.

## Reference

J.-J. E. Slotine and W. Li, *Applied Nonlinear Control*, Prentice Hall, 1991
(mass–spring–damper with nonlinear spring, p. 71).
