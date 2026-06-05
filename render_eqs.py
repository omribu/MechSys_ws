"""Render LaTeX display equations to PNG images for embedding in the .docx."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
os.makedirs("eqs", exist_ok=True)

def render(name, tex, fontsize=20):
    fig = plt.figure(figsize=(0.01,0.01))
    fig.text(0, 0, f"${tex}$", fontsize=fontsize)
    fig.savefig(f"eqs/{name}.png", dpi=200, bbox_inches='tight',
                pad_inches=0.06, transparent=True)
    plt.close(fig)

EQS = {
 # Part 1
 "mech":      r"m\,\ddot{x}(t) + b\,\dot{x}(t) + k\,x(t) + \alpha\,x^3(t) = F(t)",
 "act":       r"\tau\,\dot{F}(t) + F(t) = u(t)",
 "act_sol":   r"F(t) = u_0\left(1 - e^{-t/\tau}\right)",
 "Fexpr":     r"F(t) = m\,\ddot{x} + b\,\dot{x} + k\,x + \alpha\,x^3",
 "Fdot":      r"\dot{F}(t) = m\,\dddot{x} + b\,\ddot{x} + k\,\dot{x} + 3\alpha\,x^2\dot{x}",
 "third":     r"\tau m\,\dddot{x} + (\tau b + m)\,\ddot{x} + (\tau k + b)\,\dot{x} + k\,x + 3\tau\alpha\,x^2\dot{x} + \alpha\,x^3 = u(t)",
 "ss":        r"\dot{x}_1 = x_2,\quad \dot{x}_2 = \frac{1}{m}\left(x_3 - b x_2 - k x_1 - \alpha x_1^3\right),\quad \dot{x}_3 = \frac{1}{\tau}\left(u - x_3\right)",
 "model_form":r"y^{(n)}(t) + \sum_{i=1}^{n}[a_i\,y^{(n-i)}(t) + b_i\,f_i(y^{(n-i)}(t))] = c\,u(t)",
 "norm_div":  r"\frac{1}{c}y^{(3)} + \frac{a_1}{c}\ddot{x} + \frac{a_2}{c}\dot{x} + \frac{a_3}{c}x + \frac{b_3}{c}x^3 + \frac{g}{c}x^2\dot{x} = u",
 "norm":      r"h\,y^{(n)}(t) + \sum_{i=1}^{n}[a_i\,y^{(n-i)}(t) + b_i\,f_i(y^{(n-i)}(t))] = u(t),\qquad h=\frac{1}{c}",
 "coeffs":    r"a_1=1.25,\: a_2=1.5,\: a_3=2.0,\: b_3=1.0,\: g=1.5,\: c=2,\: h=0.5",
 # Part 2
 "refdef":    r"y_m^{(n)} + \lambda_1 y_m^{(n-1)} + \cdots + \lambda_n y_m^{(0)} = \lambda_n\, r(t)",
 "reftf":     r"y_m(s) = \dfrac{\lambda_n}{s^n + \lambda_1 s^{n-1} + \cdots + \lambda_n}\, r(s)",
 "errdef":    r"\tilde{y} = y - y_m",
 "errdyn":    r"\tilde{y}^{(n)} + \lambda_1 \tilde{y}^{(n-1)} + \cdots + \lambda_n \tilde{y}^{(0)} = 0",
 "ydes":      r"y^{(n)} = y_m^{(n)} - \lambda_1 \tilde{y}^{(n-1)} - \cdots - \lambda_n \tilde{y}^{(0)}",
 "u_known":   r"u(t) = h(y_m^{(n)} - \lambda_1\tilde{y}^{(n-1)} - \cdots - \lambda_n\tilde{y}^{(0)}) + \sum_{i=1}^{n}[a_i\,y^{(n-i)} + b_i\,f_i(y^{(n-i)})]",
 "u_adapt":   r"u(t) = h(y_m^{(n)} - \lambda_1\tilde{y}^{(n-1)} - \cdots - \lambda_n\tilde{y}^{(0)}) + \sum_{i=1}^{n}[\hat{a}_i\,y^{(n-i)} + \hat{b}_i\,f_i(y^{(n-i)})]",
 "partilde":  r"\tilde{a}_i = \hat{a}_i - a_i,\qquad \tilde{b}_i = \hat{b}_i - b_i",
 "zdef":      r"z = \tilde{y}^{(n-1)} + \lambda_1 \tilde{y}^{(n-2)} + \cdots + \lambda_{n-1}\tilde{y}^{(0)}",
 "errmodel":  r"h\,\dot{z}(t) + \lambda\,z(t) = \sum_{i=1}^{n}[\tilde{a}_i\,y^{(n-i)}(t) + \tilde{b}_i\,f_i(y^{(n-i)}(t))]",
 "adapt_law": r"\dot{\hat{a}}_i = -\gamma_i\, z(t)\, y^{(n-i)}(t),\qquad \dot{\hat{b}}_i = -\gamma_i\, z(t)\, f_i(y^{(n-i)}(t))",
 "lyap":      r"V = \frac{1}{2}(h\,z^2 + \frac{1}{\gamma}\,\tilde{\theta}^{T}\tilde{\theta})",
 "lyapdot":   r"\dot{V} = z(\tilde{\theta}^{T}\phi - \lambda z) + \frac{1}{\gamma}\,\tilde{\theta}^{T}\dot{\hat{\theta}}",
 "lyapfin":   r"\dot{V} = -\lambda\,z^2 \:\leq\: 0",
 "barbalat":  r"\ddot{V} = -2\lambda\, z\,\dot{z}\quad\text{bounded}\:\Rightarrow\:\dot{V}\to 0\:\Rightarrow\: z\to 0",
 # Part 3
 "u3_known":  r"u(t) = h(\ddot{y}_m - 2\xi\omega_n\dot{\tilde y} - \omega_n^2\tilde y - \lambda_3\tilde y \dots) + \hat{a}_1\ddot{x} + \hat{a}_2\dot{x} + \hat{a}_3 x + \hat{b}_3 x^3 + \hat{g}\,x^2\dot{x}",
 "u3":        r"u(t) = h(y_m^{(3)} - \lambda_1\tilde{y}^{(2)} - \lambda_2\tilde{y}^{(1)} - \lambda_3\tilde{y}^{(0)}) + a_1\ddot{x} + a_2\dot{x} + a_3 x + b_3 x^3 + g\,x^2\dot{x}",
 "refpoles":  r"y_m^{(3)} + \lambda_1 y_m^{(2)} + \lambda_2 y_m^{(1)} + \lambda_3 y_m = \lambda_3 r,\quad (s+p)^3,\: p=2.5",
 "reflam":    r"\lambda_1 = 3p = 7.5,\quad \lambda_2 = 3p^2 = 18.75,\quad \lambda_3 = p^3 = 15.625",
 # Part 4
 "u4":        r"u(t) = h(y_m^{(3)} - \lambda_1\tilde{y}^{(2)} - \lambda_2\tilde{y}^{(1)} - \lambda_3\tilde{y}^{(0)}) + \hat{a}_1\ddot{x} + \hat{a}_2\dot{x} + \hat{a}_3 x + \hat{b}_3 x^3 + \hat{g}\,x^2\dot{x}",
 "adapt4":    r"\dot{\hat{a}}_1=-\gamma z\,\ddot{x},\:\: \dot{\hat{a}}_2=-\gamma z\,\dot{x},\:\: \dot{\hat{a}}_3=-\gamma z\,x,\:\: \dot{\hat{b}}_3=-\gamma z\,x^3,\:\: \dot{\hat{g}}=-\gamma z\,x^2\dot{x}",
 "zdef3":     r"z = \tilde{y}^{(2)} + \lambda_1 \tilde{y}^{(1)} + \lambda_2 \tilde{y}^{(0)}",
 # Part 5
 "bdiff":     r"\dot{y}(k) \approx \dfrac{y(k) - y(k-1)}{T}",
 "dirty":     r"\dot{y}_f(k) = a_f\,\dot{y}_f(k-1) + (1-a_f)\,\dfrac{y(k)-y(k-1)}{T},\qquad a_f = \dfrac{\tau_d}{\tau_d + T}",
 "zoh":       r"u(t) = u(kT),\qquad kT \leq t < (k+1)T",
 "adapt_disc":r"\hat{\theta}(k) = \hat{\theta}(k-1) - \gamma\,\dfrac{z(k)\,\phi(k)}{1+\phi^{T}\phi}\,T",
 # Part 6
 "utotal":    r"u_{\text{total}} = u_{\text{mrac}} + u_{\text{pid}}",
 "upid":      r"u_{\text{pid}} = K_P\,\tilde{y} + K_I\int_0^t \tilde{y}\,d\tau + K_D\,\dot{\tilde{y}}",
 "pidgains":  r"K_P = 0.8,\quad K_I = 0.1,\quad K_D = 6.0",
}

for name, tex in EQS.items():
    render(name, tex)
print(f"Rendered {len(EQS)} equation images.")
