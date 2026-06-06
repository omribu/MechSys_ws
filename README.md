# Mechatronic Systems HW —  LaTeX 

Nonlinear adaptive control (MRAC) of a third-order mass–spring–damper with a
nonlinear cubic spring and a first-order actuator. 

## One convention everywhere

- Parameter vector  θ = [a₁, a₂, a₃, b₁, b₂]
- Regressor         φ = [ÿ, ẏ, y, y²ẏ, y³]
- Nonlinear funcs   f₁ = y²ẏ (mixed term),  f₂ = y³ (cubic term)
  → so **b₁ = 3α/m = 3.0** (mixed),  **b₂ = α/(τm) = 2.0** (cubic)
- Plant in standard form:  y⁽³⁾ + a₁ÿ + a₂ẏ + a₃y + b₁f₁ + b₂f₂ = c·u   (c = 2.0)
- Filtered error    s = ë + α₂ė + α₁e        (e = y − y_m)
- Controller gains  k₂ = α₂+λ, k₁ = α₁+λα₂, k₀ = λα₁
- Reference model   y_m⁽³⁾ + λ₂ÿ_m + λ₁ẏ_m + λ₀y_m = λ₀r   (λ₀,λ₁,λ₂ Hurwitz)
- Adaptation        θ̂̇ = −Γ φ s


## Files

| file/folder          | description                                      |
|----------------------|--------------------------------------------------|
| `main.tex`           | The full report (LaTeX), Sections 1–6.           |
| `main.pdf`           | Pre-compiled PDF (22 pages).                      |
| `mechatronics_hw.py` | Single Python file that generates every figure.  |
| `figures/`, `figs/`  | All figures (both paths are on `\graphicspath`). |


## Regenerate figures

```
pip install numpy matplotlib
python3 mechatronics_hw.py     # writes ./figures/
```


## Numbers in the text match the figures (simulation)

- θ = [2.5, 3.0, 4.0, 3.0, 2.0], c = 2.0
- Reference: triple pole p = 2.5 → λ₂=7.5, λ₁=18.75, λ₀=15.625
- Gains: α₁=6, α₂=5, λ=5 → k₂=10, k₁=31, k₀=30
- Nominal tracking error ≈ 0 (machine precision); adaptive ≈ 3×10⁻⁷
- Continuous RMSE 0.035; discrete 20 Hz RMSE 0.048
- Sampling sweep: 0.048 (20 ms) → 0.053 (65 ms), unstable at 70 ms
- PID at 50 ms: 0.047 vs 0.048

## Reference

J.-J. E. Slotine and W. Li, *Applied Nonlinear Control*, Prentice Hall, 1991
(mass–spring–damper with nonlinear spring, p. 71).
