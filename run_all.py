"""
run_all.py — Master driver for the Mechatronic Systems homework.

Runs every simulation and regenerates all figures into ./figs/.
Just execute:  python3 run_all.py

Covers:
  Part 3 : nominal control with known parameters        (sim_p34.py)
  Part 4 : adaptive control with unknown parameters      (sim_p34.py)
  Part 5 : indirect digital (discrete-time) controller   (sim_p56.py)
  Part 6 : sampling-time degradation + PID remedy        (sim_p56.py)
  Diagrams: physical schematic, MRAC block diagram, etc. (make_diagrams.py)
"""
import os, runpy

os.makedirs("figs", exist_ok=True)

print("="*60)
print("Mechatronic Systems HW — running all simulations")
print("="*60)

for name in ["make_diagrams.py", "sim_p34.py", "sim_p56.py"]:
    print(f"\n>>> running {name}")
    runpy.run_path(name, run_name="__main__")

print("\n" + "="*60)
print("Done. All figures are in ./figs/")
print("="*60)
