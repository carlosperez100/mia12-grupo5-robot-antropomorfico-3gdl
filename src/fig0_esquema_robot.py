"""
Figura 0. Esquema del brazo RRR del paper (vista lateral en el plano del brazo y vista superior),
con las longitudes de la Tabla 3 y los ángulos θ1, θ2 y θ3.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from pathlib import Path
from modelo_robot import L1, L2, L3

FIG = Path(__file__).resolve().parents[1] / "assets" / "figuras"
AZUL, ROJO, GRIS = "#1F3A5F", "#C0392B", "#5B6B7C"

t2, t3 = np.radians(35), np.radians(-70)
hombro = np.array([0.0, L1])
codo = hombro + L2 * np.array([np.cos(t2), np.sin(t2)])
punta = codo + L3 * np.array([np.cos(t2 + t3), np.sin(t2 + t3)])

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.5, 1]})

# Vista lateral
ax.fill_between([-0.06, 0.06], -0.02, 0, color=GRIS, alpha=0.4)
ax.plot([0, 0], [0, L1], "-", color=AZUL, lw=8, solid_capstyle="round")
ax.plot([hombro[0], codo[0]], [hombro[1], codo[1]], "-", color=AZUL, lw=8, solid_capstyle="round")
ax.plot([codo[0], punta[0]], [codo[1], punta[1]], "-", color=AZUL, lw=8, solid_capstyle="round")
for p in (hombro, codo):
    ax.plot(*p, "o", ms=14, color="white", mec=AZUL, mew=3)
ax.plot(*punta, "o", ms=10, color=ROJO)
ax.text(0.25, 0.265, "punta P = (x, y, z)", color=ROJO, fontsize=10)
ax.text(-0.03, L1 / 2, "l1 = 0.32 m\n(base)", fontsize=10, ha="right", va="center")
m2 = (hombro + codo) / 2
ax.text(m2[0] - 0.02, m2[1] + 0.03, "l2 = 0.22 m\n(brazo)", fontsize=10, ha="right")
m3 = (codo + punta) / 2
ax.text(0.30, 0.435, "l3 = 0.22 m\n(antebrazo)", fontsize=10, va="top")
ax.plot([hombro[0], hombro[0] + 0.12], [hombro[1], hombro[1]], "--", color=GRIS, lw=1)
ax.add_patch(Arc(hombro, 0.16, 0.16, theta1=0, theta2=np.degrees(t2), color=ROJO, lw=2))
ax.text(hombro[0] + 0.09, hombro[1] + 0.03, "θ2 (hombro)", color=ROJO, fontsize=10)
ext = codo + 0.1 * np.array([np.cos(t2), np.sin(t2)])
ax.plot([codo[0], ext[0]], [codo[1], ext[1]], "--", color=GRIS, lw=1)
ax.add_patch(Arc(codo, 0.14, 0.14, theta1=np.degrees(t2 + t3), theta2=np.degrees(t2), color=ROJO, lw=2))
ax.text(codo[0] + 0.08, codo[1] + 0.03, "θ3 (codo)", color=ROJO, fontsize=10)
ax.annotate("", xy=(0.0, 0.05), xytext=(0.0, -0.0), arrowprops=dict(arrowstyle="-", color=GRIS))
ax.set_xlim(-0.25, 0.5); ax.set_ylim(-0.03, 0.62); ax.set_aspect("equal")
ax.set_xlabel("r = distancia horizontal desde el eje de la base [m]"); ax.set_ylabel("z [m]")
ax.set_title("Vista lateral (plano del brazo)", color=AZUL)
ax.grid(alpha=0.25)

# Vista superior
t1 = np.radians(30)
r = punta[0]
ax2.plot([0, r * np.cos(t1)], [0, r * np.sin(t1)], "-", color=AZUL, lw=6, solid_capstyle="round")
ax2.plot(0, 0, "o", ms=14, color="white", mec=AZUL, mew=3)
ax2.plot(r * np.cos(t1), r * np.sin(t1), "o", ms=10, color=ROJO)
ax2.plot([0, 0.3], [0, 0], "--", color=GRIS, lw=1)
ax2.add_patch(Arc((0, 0), 0.2, 0.2, theta1=0, theta2=30, color=ROJO, lw=2))
ax2.text(0.11, 0.025, "θ1 (base)", color=ROJO, fontsize=10)
ax2.text(0.24, -0.025, "eje x0", fontsize=9, color=GRIS)
ax2.text(-0.09, -0.075, "la base gira alrededor del eje z0 (vertical)", fontsize=9, color=GRIS)
ax2.set_xlim(-0.1, 0.35); ax2.set_ylim(-0.1, 0.25); ax2.set_aspect("equal"); ax2.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
ax2.set_xlabel("x [m]"); ax2.set_ylabel("y [m]")
ax2.set_title("Vista superior", color=AZUL)
ax2.grid(alpha=0.25)

fig.suptitle("Brazo antropomórfico RRR del paper: base, hombro y codo (Tablas 2 y 3)", color=AZUL, fontsize=13)
fig.tight_layout()
fig.savefig(FIG / "fig0_esquema_robot.png", dpi=150)
print("ok")
