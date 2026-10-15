"""
Paso 3. Espacio de trabajo con las longitudes (Tabla 3) y límites articulares (Tabla 1) del paper.
Mismo procedimiento de clase: barrido de las articulaciones + cinemática directa.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from modelo_robot import L1, L2, L3, LIMITES

ROOT = Path(__file__).resolve().parents[1]
OUT, FIG = ROOT / "resultados", ROOT / "assets" / "figuras"

n = 45
t1 = np.radians(np.linspace(*LIMITES["theta1"], n))
t2 = np.radians(np.linspace(*LIMITES["theta2"], n))
t3 = np.radians(np.linspace(*LIMITES["theta3"], n))
T1, T2, T3 = np.meshgrid(t1, t2, t3, indexing="ij")
R = L2 * np.cos(T2) + L3 * np.cos(T2 + T3)
X, Y, Z = np.cos(T1) * R, np.sin(T1) * R, L1 + L2 * np.sin(T2) + L3 * np.sin(T2 + T3)

# Corte en el plano del brazo (θ1 = 0): r vs z
rng = np.random.default_rng(1)
t2f = np.radians(rng.uniform(*LIMITES["theta2"], 60000))
t3f = np.radians(rng.uniform(*LIMITES["theta3"], 60000))
rf = L2 * np.cos(t2f) + L3 * np.cos(t2f + t3f)
zf = L1 + L2 * np.sin(t2f) + L3 * np.sin(t2f + t3f)

seqA = np.array([(0.25, 0.01, 0.11), (0.25, 0.11, 0.11), (0.15, 0.11, 0.11), (0.15, 0.01, 0.11)])

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(1, 2, 1, projection="3d")
idx = np.random.default_rng(0).choice(X.size, 25000, replace=False)
ax.scatter(X.ravel()[idx], Y.ravel()[idx], Z.ravel()[idx], s=0.4, c=Z.ravel()[idx], cmap="viridis", alpha=0.5)
ax.scatter(seqA[:, 0], seqA[:, 1], seqA[:, 2], c="red", s=30, label="Puntos de la Tabla 10")
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
ax.set_title("Espacio de trabajo 3D (límites de la Tabla 1)")
ax.legend(fontsize=8); ax.set_box_aspect((1, 1, 0.8))

ax2 = fig.add_subplot(1, 2, 2)
ax2.scatter(rf, zf, s=0.4, c="#2F6FB0", alpha=0.25)
ax2.scatter(np.hypot(seqA[:, 0], seqA[:, 1]), seqA[:, 2], c="red", s=30, label="Puntos de la Tabla 10")
ax2.plot(0, L1, "ks", label="Hombro (z = l1)")
ax2.set_xlabel("r = √(x² + y²) [m]"); ax2.set_ylabel("z [m]")
ax2.set_title("Corte en el plano del brazo")
ax2.set_aspect("equal"); ax2.grid(alpha=0.3); ax2.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "fig3_espacio_trabajo.png", dpi=150)

res = {
    "alcance_max_desde_hombro_m": L2 + L3,
    "r_max_m": round(float(R.max()), 4),
    "z_min_m": round(float(Z.min()), 4), "z_max_m": round(float(Z.max()), 4),
    "puntos_tabla10_distancia_al_hombro_m": [round(float(np.hypot(np.hypot(p[0], p[1]), p[2] - L1)), 4) for p in seqA],
}
(OUT / "paso3_espacio_trabajo.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
print(res)
