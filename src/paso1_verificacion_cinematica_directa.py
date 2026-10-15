"""
Paso 1. Comprobación de la cinemática directa publicada en el paper.

Se compara el producto de matrices DH (Tabla 2 del paper, armado con el
método del curso) contra las ecuaciones (9), (12), (13) y (14) tal como
aparecen impresas en el artículo.
"""
import json
import numpy as np
import sympy as sp
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "resultados"
OUT.mkdir(exist_ok=True)

t1, t2, t3, l1, l2, l3 = sp.symbols("theta1 theta2 theta3 l1 l2 l3", real=True)
C1, S1, C2, S2, C3, S3 = sp.cos(t1), sp.sin(t1), sp.cos(t2), sp.sin(t2), sp.cos(t3), sp.sin(t3)


def dh(theta, d, a, alpha):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta) * sp.cos(alpha),  sp.sin(theta) * sp.sin(alpha), a * sp.cos(theta)],
        [sp.sin(theta),  sp.cos(theta) * sp.cos(alpha), -sp.cos(theta) * sp.sin(alpha), a * sp.sin(theta)],
        [0, sp.sin(alpha), sp.cos(alpha), d],
        [0, 0, 0, 1],
    ])


A1 = dh(t1, l1, 0, sp.pi / 2)
A2 = dh(t2, 0, l2, 0)
A3 = dh(t3, 0, l3, 0)
T = sp.simplify(A1 * A2 * A3)
P = T[:3, 3]

# --- Ecuación (9) del paper: 1A2 impresa con "l2C2" en la fila 2, columna 4
A2_paper = sp.Matrix([[C2, -S2, 0, l2 * C2],
                      [S2,  C2, 0, l2 * C2],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
dif_eq9 = sp.simplify(A2_paper - A2)

# --- Ecuaciones (12)-(14) del paper, transcritas tal cual
Px_paper = l3 * C1 * C2 * C3 - l3 * C1 * S2 * S3 + l2 * C1 * C2
Py_paper = l3 * S3 * C2 * C3 - l3 * S1 * S2 * S3 + l3 * S1 * C2
Pz_paper = l3 * S2 * C3 + l3 * C2 * S3 + l2 * S2 + l1

res = {
    "posicion_DH_simplificada": [str(sp.simplify(e)) for e in P],
    "ec9_diferencia_paper_menos_DH": str(dif_eq9.tolist()),
    "ec12_Px_coincide": sp.simplify(Px_paper - P[0]) == 0,
    "ec13_Py_coincide": sp.simplify(Py_paper - P[1]) == 0,
    "ec13_Py_diferencia": str(sp.simplify(Py_paper - P[1])),
    "ec14_Pz_coincide": sp.simplify(Pz_paper - P[2]) == 0,
}

# Ejemplo numérico para mostrar el efecto de la errata de la ec. (13)
vals = {l1: 0.32, l2: 0.22, l3: 0.22, t1: sp.rad(30), t2: sp.rad(20), t3: sp.rad(-60)}
res["ejemplo_q_grados"] = [30, 20, -60]
res["ejemplo_P_DH"] = [round(float(e.subs(vals)), 4) for e in P]
res["ejemplo_Py_paper"] = round(float(Py_paper.subs(vals)), 4)

(OUT / "paso1_cinematica_directa.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
for k, v in res.items():
    print(f"{k}: {v}")

# --- Figura: Py por DH vs Py con la ec. (13) impresa, barriendo θ1 con θ2 = 20°, θ3 = -60°
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
f_dh = sp.lambdify(t1, P[1].subs({l1: 0.32, l2: 0.22, l3: 0.22, t2: sp.rad(20), t3: sp.rad(-60)}), "numpy")
f_pa = sp.lambdify(t1, Py_paper.subs({l1: 0.32, l2: 0.22, l3: 0.22, t2: sp.rad(20), t3: sp.rad(-60)}), "numpy")
th = np.linspace(-155, 155, 300)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(th, f_dh(np.radians(th)), lw=2.5, color="#1F3A5F", label="Py por producto DH (Tabla 2)")
ax.plot(th, f_pa(np.radians(th)) * np.ones_like(th), "--", lw=2, color="#C0392B", label="Py con la ec. (13) impresa")
ax.set_xlabel("θ1 [°]  (θ2 = 20°, θ3 = −60°)"); ax.set_ylabel("Py [m]")
ax.set_title("Errata de la ec. (13): la coordenada y del extremo")
ax.grid(alpha=0.3); ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(Path(__file__).resolve().parents[1] / "assets" / "figuras" / "fig1_errata_ec13.png", dpi=150)
