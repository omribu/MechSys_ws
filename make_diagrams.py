"""Generate the physical-system schematic, free-body diagram, and the MRAC block
diagram used in Parts 1 and 2 of the report."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import numpy as np

plt.rcParams.update({"font.size": 12})

# ---------------- Physical system + FBD ----------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.6))

def wall(ax, x, y0, y1):
    ax.plot([x, x], [y0, y1], 'k', lw=2)
    for yy in np.linspace(y0, y1, 9):
        ax.plot([x-0.18, x], [yy, yy+0.12], 'k', lw=1)

def spring(ax, x0, x1, y, n=7, amp=0.12):
    xs = np.linspace(x0, x1, 2*n+1)
    ys = y + amp*np.array([0]+[(-1)**i for i in range(2*n-1)]+[0])
    ax.plot(xs, ys, 'k', lw=1.4)

def nlspring(ax, x0, x1, y, n=9, amp=0.16):
    xs = np.linspace(x0, x1, 200)
    ys = y + amp*np.sin(np.linspace(0, n*np.pi, 200))*np.hanning(200)*1.6
    ax.plot(xs, ys, 'k', lw=1.4)

def damper(ax, x0, x1, y):
    xm=(x0+x1)/2
    ax.plot([x0, xm-0.12], [y, y], 'k', lw=1.4)
    ax.add_patch(mp.Rectangle((xm-0.12, y-0.16), 0.24, 0.32, fill=False, lw=1.4))
    ax.plot([xm, x1], [y, y], 'k', lw=1.4)
    ax.plot([xm-0.06, xm-0.06], [y-0.12, y+0.12], 'k', lw=2)

# Left: physical system
axL.set_title("Physical system", fontsize=13, fontweight='bold')
wall(axL, 0.0, 0.2, 3.0)
yk, yb, yn = 2.4, 1.6, 0.8
spring(axL, 0.0, 2.0, yk); axL.text(1.0, yk+0.3, "$k$", fontsize=14)
damper(axL, 0.0, 2.0, yb); axL.text(1.0, yb+0.3, "$b$", fontsize=14)
nlspring(axL, 0.0, 2.0, yn); axL.text(0.55, yn-0.55, "nonlinear spring: $\\alpha x^3$", fontsize=10)
mass = mp.Rectangle((2.0, 0.45), 0.9, 2.5, facecolor="#e8e8e8", edgecolor='k', lw=1.6)
axL.add_patch(mass); axL.text(2.45, 1.7, "$m$", fontsize=16, ha='center')
# wheels
for wx in [2.2, 2.7]:
    axL.add_patch(mp.Circle((wx, 0.33), 0.12, fill=False, lw=1.4))
axL.plot([1.9, 3.4], [0.18, 0.18], 'k', lw=1.5)
for xx in np.linspace(1.95, 3.35, 12):
    axL.plot([xx, xx-0.1], [0.18, 0.06], 'k', lw=0.8)
# actuator + force
axL.add_patch(mp.Rectangle((2.9, 1.45), 0.28, 0.3, facecolor="#cccccc", edgecolor='k'))
axL.annotate("", xy=(4.3, 1.6), xytext=(3.2, 1.6),
             arrowprops=dict(arrowstyle="-|>", color="C0", lw=2.2))
axL.text(3.95, 1.85, "$F(t)$", color="C0", fontsize=14)
axL.text(3.0, 1.15, "actuator", color="C0", fontsize=10)
axL.annotate("", xy=(3.3, 0.0), xytext=(2.45, 0.0),
             arrowprops=dict(arrowstyle="-|>", color="C0", lw=1.8))
axL.text(2.8, -0.28, "$x(t)$", color="C0", fontsize=13)
axL.set_xlim(-0.4, 4.6); axL.set_ylim(-0.6, 3.4); axL.axis('off')

# Right: free body diagram
axR.set_title("Free-body diagram of the mass", fontsize=13, fontweight='bold')
mass2 = mp.Rectangle((2.0, 0.45), 0.9, 2.5, facecolor="#e8e8e8", edgecolor='k', lw=1.6)
axR.add_patch(mass2); axR.text(2.45, 1.7, "$m$", fontsize=16, ha='center')
forces = [("$-kx$", 2.45), ("$-b\\dot x$", 1.85), ("$-\\alpha x^3$", 1.25)]
for lbl, yy in forces:
    axR.annotate("", xy=(0.9, yy), xytext=(2.0, yy),
                 arrowprops=dict(arrowstyle="-|>", color='k', lw=1.8))
    axR.text(0.55, yy+0.02, lbl, fontsize=13, ha='right')
axR.annotate("", xy=(4.2, 1.7), xytext=(2.9, 1.7),
             arrowprops=dict(arrowstyle="-|>", color="C0", lw=2.2))
axR.text(3.85, 1.95, "$F(t)$", color="C0", fontsize=14)
axR.annotate("", xy=(3.3, 0.0), xytext=(2.45, 0.0),
             arrowprops=dict(arrowstyle="-|>", color="C0", lw=1.8))
axR.text(2.8, -0.28, "$x(t)$", color="C0", fontsize=13)
axR.set_xlim(-0.4, 4.6); axR.set_ylim(-0.6, 3.4); axR.axis('off')
plt.tight_layout(); plt.savefig("figs/schematic.png", dpi=130); plt.close()

# ---------------- MRAC block diagram ----------------
fig, ax = plt.subplots(figsize=(10, 4.4))
def box(x, y, w, h, text, fc="#eef3fb"):
    ax.add_patch(mp.FancyBboxBox if False else mp.Rectangle((x,y),w,h, facecolor=fc, edgecolor='k', lw=1.5))
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=11)
def arrow(x0,y0,x1,y1,txt=None,off=0.12,color='k'):
    ax.annotate("", xy=(x1,y1), xytext=(x0,y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6))
    if txt: ax.text((x0+x1)/2, (y0+y1)/2+off, txt, ha='center', fontsize=10)

# reference model
box(0.3, 2.7, 1.7, 0.8, "Reference\nmodel $W_m(s)$", fc="#e7f0e7")
# controller
box(3.2, 2.7, 2.0, 0.8, "Adaptive\ncontroller", fc="#eef3fb")
# plant
box(6.6, 2.7, 1.9, 0.8, "Nonlinear\nplant (3rd order)", fc="#fbeeee")
# adaptation
box(3.2, 0.7, 2.0, 0.8, "Adaptation law\n$\\dot{\\hat\\theta}=-\\gamma z\\,\\phi$", fc="#f5f0e0")

arrow(-0.4, 3.1, 0.3, 3.1, "$r(t)$")
arrow(2.0, 3.1, 3.2, 3.1, "$y_m$")
arrow(5.2, 3.1, 6.6, 3.1, "$u(t)$")
arrow(8.5, 3.1, 9.4, 3.1, "$y(t)$")
# feedback of y to controller
ax.plot([9.1,9.1],[3.1,1.9],'k',lw=1.4); ax.plot([9.1,2.6],[1.9,1.9],'k',lw=1.4)
ax.annotate("",xy=(2.6,2.7),xytext=(2.6,1.9),arrowprops=dict(arrowstyle="-|>",color='k',lw=1.4))
# error to adaptation
ax.plot([1.15,1.15],[2.7,1.1],'k',lw=1.4)
ax.annotate("",xy=(3.2,1.1),xytext=(1.15,1.1),arrowprops=dict(arrowstyle="-|>",color='k',lw=1.4))
ax.text(1.6,1.2,"$\\tilde y = y-y_m$",fontsize=10)
# y also into error node (summing): show y feeding the e line
ax.plot([9.1,9.1],[1.9,0.5],'k',lw=1.0,ls=':')
ax.plot([9.1,1.15],[0.5,0.5],'k',lw=1.0,ls=':')
ax.plot([1.15,1.15],[0.5,1.1],'k',lw=1.0,ls=':')
# adaptation updates controller
ax.annotate("",xy=(4.2,2.7),xytext=(4.2,1.5),arrowprops=dict(arrowstyle="-|>",color='C3',lw=1.6))
ax.text(4.3,2.0,"$\\hat\\theta$",color='C3',fontsize=12)
ax.set_xlim(-0.6,9.6); ax.set_ylim(0.3,3.8); ax.axis('off')
ax.set_title("Model-Reference Adaptive Control (MRAC) structure", fontsize=12, fontweight='bold')
plt.tight_layout(); plt.savefig("figs/blockdiagram.png", dpi=130); plt.close()

# ---------------- Reference-model step response ----------------
from sim_core import L1,L2,L3
import numpy as np
dt=1e-3; t=np.arange(0,5,dt); ym=np.zeros(len(t)); ymd=0; ymdd=0; ymv=0
Y=np.zeros(len(t))
for i in range(len(t)):
    ymddd=L3*1.0-L1*ymdd-L2*ymd-L3*ymv
    Y[i]=ymv
    ymv+=ymd*dt; ymd+=ymdd*dt; ymdd+=ymddd*dt
plt.figure(figsize=(6.5,3.6))
plt.plot(t,Y,'C0',lw=1.8); plt.axhline(1.0,color='k',ls='--',lw=0.8)
plt.xlabel('t [s]'); plt.ylabel('$y_m$'); plt.title('Reference-model step response (triple pole at $p=2.5$)')
plt.grid(alpha=0.3); plt.tight_layout(); plt.savefig('figs/refmodel.png',dpi=130); plt.close()

print("Schematic, block diagram, reference-model figures written.")
