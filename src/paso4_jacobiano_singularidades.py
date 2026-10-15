"""
Paso 4. Jacobiano y singularidades (no desarrollados en el paper).

- Jacobiano analítico: derivada simbólica de la posición respecto a (θ1, θ2, θ3).
- Jacobiano geométrico: Jv_i = z_{i-1} × (p - p_{i-1}), Jw_i = z_{i-1}. Ambos métodos vistos en clase.
- Singularidades: det(Jv) = 0.
- Evaluación numérica en los puntos de la Tabla 10.
"""
import json
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from modelo_robot import cinematica_inversa, jacobiano, L1, L2, L3

ROOT = Path(__file__).resolve().parents[1]
OUT, FIG = ROOT / "resultados", ROOT / "assets" / "figuras"

t1, t2, t3, l1, l2, l3 = sp.symbols("theta1 theta2 theta3 l1 l2 l3", real=True)


def dh(theta, d, a, alpha):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta) * sp.cos(alpha),  sp.sin(theta) * sp.sin(alpha), a * sp.cos(theta)],
        [sp.sin(theta),  sp.cos(theta) * sp.cos(alpha), -sp.cos(theta) * sp.sin(alpha), a * sp.sin(theta)],
        [0, sp.sin(alpha), sp.cos(alpha), d], [0, 0, 0, 1]])


A1, A2, A3 = dh(t1, l1, 0, sp.pi / 2), dh(t2, 0, l2, 0), dh(t3, 0, l3, 0)
T01, T02, T03 = A1, sp.simplify(A1 * A2), sp.simplify(A1 * A2 * A3)
P = T03[:3, 3]
q = sp.Matrix([t1, t2, t3])

# Jacobiano analítico (derivada directa)
J_an = sp.simplify(P.jacobian(q))

# Jacobiano geométrico
z0, p0 = sp.Matrix([0, 0, 1]), sp.Matrix([0, 0, 0])
z1, p1 = T01[:3, 2], T01[:3, 3]
z2, p2 = T02[:3, 2], T02[:3, 3]
Jv = sp.simplify(sp.Matrix.hstack(z0.cross(P - p0), z1.cross(P - p1), z2.cross(P - p2)))
Jw = sp.Matrix.hstack(z0, z1, z2)

coinciden = sp.simplify(J_an - Jv) == sp.zeros(3, 3)
det = sp.factor(sp.simplify(J_an.det()))
print("Jacobiano analítico = geométrico (parte lineal):", coinciden)
print("det(Jv) =", det)
print("Jw =", Jw.T.tolist())

# Evaluación en los puntos de la Tabla 10 (solución codo arriba)
seqA = [(0.25, 0.01, 0.11), (0.25, 0.11, 0.11), (0.15, 0.11, 0.11), (0.15, 0.01, 0.11)]
evals = []
for p in seqA:
    qq = cinematica_inversa(p)[0]
    J = jacobiano(qq)
    sv = np.linalg.svd(J, compute_uv=False)
    evals.append({"punto": p, "det": round(float(np.linalg.det(J)), 6),
                  "numero_condicion": round(float(sv[0] / sv[-1]), 2)})
for e in evals:
    print(e)

# Mapa de det(J) en el plano (θ2, θ3)
g2, g3 = np.meshgrid(np.linspace(-35, 130, 300), np.linspace(-130, 130, 300), indexing="ij")
D = -L2 * L3 * np.sin(np.radians(g3)) * (L2 * np.cos(np.radians(g2)) + L3 * np.cos(np.radians(g2 + g3)))
fig, ax = plt.subplots(figsize=(6.5, 5))
cs = ax.contourf(g2, g3, D, levels=30, cmap="RdBu_r")
ax.contour(g2, g3, D, levels=[0], colors="k", linewidths=2)
fig.colorbar(cs, label="det(J) [m³]")
for p in seqA:
    qq = np.degrees(cinematica_inversa(p)[0])
    ax.plot(qq[1], qq[2], "o", color="gold", mec="k")
ax.set_xlabel("θ2 [°]"); ax.set_ylabel("θ3 [°]")
ax.set_title("det(J) y singularidades (línea negra = det(J) = 0)")
fig.tight_layout()
fig.savefig(FIG / "fig4_singularidades.png", dpi=150)

res = {"jacobiano_analitico": str(J_an.tolist()), "Jw": str(Jw.tolist()),
       "analitico_igual_geometrico": bool(coinciden), "det_Jv": str(det),
       "evaluacion_tabla10": evals}
(OUT / "paso4_jacobiano.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
