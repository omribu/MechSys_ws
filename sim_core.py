"""
Mechatronic Systems HW - 3rd order nonlinear mass-spring-damper + 1st order actuator.

Plant (third order, after eliminating actuator force F):
    tau*m x''' + (tau*b+m) x'' + (tau*k+b) x' + k x + 3*tau*alpha x^2 x' + alpha x^3 = u

Divided by (tau*m) -> standard form  y''' + a1 y'' + a2 y' + a3 y + b3 y^3 + g*y^2*y' = c u
with y = x.

State-space (x1=x, x2=x', x3=F):
    x1' = x2
    x2' = (1/m)(x3 - b x2 - k x1 - alpha x1^3)
    x3' = (1/tau)(u - x3)
"""
import numpy as np

# ----------------------- physical parameters -----------------------
M, B, K, ALPHA, TAU = 1.0, 0.5, 2.0, 1.0, 0.5

# Physical (un-normalized) model-form coefficients, y = x, n = 3:
#   y''' + A1p y'' + A2p y' + A3p y + B3p y^3 + Gp y^2 y' = C u
A1p = (TAU*B + M)/(TAU*M)     # coeff of x''   = 2.5
A2p = (TAU*K + B)/(TAU*M)     # coeff of x'    = 3.0
A3p = K/(TAU*M)              # coeff of x     = 4.0
B3p = ALPHA/(TAU*M)          # coeff of x^3   = 2.0
Gp  = 3*ALPHA/M              # coeff of x^2 x'= 3.0
C   = 1.0/(TAU*M)            # input gain c   = 2.0

# Reference-style normalization: divide the whole equation by c so that
#   h y''' + a1 y'' + a2 y' + a3 y + b3 y^3 + g y^2 y' = u,  with h = 1/c.
H  = 1.0/C                   # 0.5
A1 = A1p/C                   # 1.25
A2 = A2p/C                   # 1.5
A3 = A3p/C                   # 2.0
B3 = B3p/C                   # 1.0
G  = Gp/C                    # 1.5


def plant_step(state, u, dt):
    """One Euler step of the true third-order plant in physical state-space form."""
    x1, x2, x3 = state
    dx1 = x2
    dx2 = (1.0/M)*(x3 - B*x2 - K*x1 - ALPHA*x1**3)
    dx3 = (1.0/TAU)*(u - x3)
    return np.array([x1 + dx1*dt, x2 + dx2*dt, x3 + dx3*dt])


def regressors(y, yd, ydd):
    """
    Known regressor functions phi(state) that multiply the (a_i, b_i) gains
    in normalized model form. y=x, yd=x', ydd=x''.
    Order: [x'' , x' , x , x^3 , x^2 x']  matching [a1, a2, a3, b3, g].
    """
    return np.array([ydd, yd, y, y**3, (y**2)*yd])


# Reference model: 3rd order, stable, Hurwitz.  desired closed-loop poles.
# y_m''' + L1 y_m'' + L2 y_m' + L3 y_m = L3 r   (unity DC gain)
# choose triple real pole at p -> (s+p)^3
P_REF = 2.5
L1 = 3*P_REF
L2 = 3*P_REF**2
L3 = P_REF**3
LAM = np.array([L1, L2, L3])   # error-dynamics gains lambda_1..lambda_3


def ref_step(ym_state, r, dt):
    """Euler step of reference model. ym_state = [ym, ym', ym'']. returns new state + ym'''."""
    ym, ymd, ymdd = ym_state
    ymddd = L3*r - L1*ymdd - L2*ymd - L3*ym
    return np.array([ym + ymd*dt, ymd + ymdd*dt, ymdd + ymddd*dt]), ymddd
