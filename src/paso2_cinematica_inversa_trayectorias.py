"""
Paso 2. Cinemática inversa geométrica sobre los puntos publicados en el paper.

Tablas 10 y 11 del paper: coordenadas cartesianas (m) de las secuencias A y B.
Para cada punto se calculan las dos soluciones (codo arriba / codo abajo), se
comprueba con la cinemática directa y se revisan los límites articulares.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from modelo_robot import (cinematica_inversa, cinematica_directa, dentro_de_limites,
                          det_jacobiano, matrices_dh, L1)

ROOT = Path(__file__).resolve().parents[1]
OUT, FIG = ROOT / "resultados", ROOT / "assets" / "figuras"

SECUENCIA_A = [(0.25, 0.01, 0.11), (0.25, 0.01, 0.11), (0.25, 0.11, 0.11),
               (0.15, 0.11, 0.11), (0.15, 0.01, 0.11), (0.25, 0.01, 0.11)]
SECUENCIA_B = [(0.25, 0.01, 0.11), (0.25, 0.01, 0.11), (0.25, 0.11, 0.11),
               (0.15, 0.11, 0.11), (0.25, 0.01, 0.11), (0.25, 0.01, 0.11)]

filas = []
for nombre, seq in (("A", SECUENCIA_A), ("B", SECUENCIA_B)):
    for k, p in enumerate(seq, start=1):
        for etiqueta, q in zip(("codo arriba", "codo abajo"), cinematica_inversa(p)):
            p_fk = cinematica_directa(q)
            err = np.linalg.norm(p_fk - np.array(p))
            g = np.degrees(q)
            filas.append({
                "secuencia": nombre, "punto": k, "x": p[0], "y": p[1], "z": p[2],
                "solucion": etiqueta,
                "theta1_deg": round(g[0], 3), "theta2_deg": round(g[1], 3), "theta3_deg": round(g[2], 3),
                "error_FK_m": f"{err:.2e}",
                "dentro_limites": dentro_de_limites(q),
                "det_J": round(det_jacobiano(q), 6),
            })

with open(OUT / "paso2_cinematica_inversa_tablas10_11.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=filas[0].keys())
    w.writeheader(); w.writerows(filas)

for r in filas:
    if r["secuencia"] == "A":
        print(r)

# --- Figura: robot en cada punto de la secuencia A (codo arriba)
def puntos_brazo(q):
    A1, A2, A3 = matrices_dh(q)
    T01, T02, T03 = A1, A1 @ A2, A1 @ A2 @ A3
    return np.array([[0, 0, 0], T01[:3, 3], T02[:3, 3], T03[:3, 3]])

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")
seq = np.array(SECUENCIA_A)
ax.plot(seq[:, 0], seq[:, 1], seq[:, 2], "k--", lw=1.5, label="Trayectoria A (Tabla 10 del paper)")
colores = plt.cm.viridis(np.linspace(0, 1, len(SECUENCIA_A)))
for k, (p, c) in enumerate(zip(SECUENCIA_A, colores), start=1):
    q = cinematica_inversa(p)[0]
    b = puntos_brazo(q)
    ax.plot(b[:, 0], b[:, 1], b[:, 2], "-o", color=c, lw=2.5, ms=4)
    if k in (1, 3, 4, 5):
        etiqueta = "P1 = P2 = P6" if k == 1 else f"P{k}"
        ax.text(p[0], p[1], p[2] - 0.03, etiqueta, fontsize=8)
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
ax.set_title("Cinemática inversa en los puntos de la secuencia A (codo arriba)")
ax.legend(loc="upper left", fontsize=8)
ax.set_box_aspect((1, 1, 1))
ax.set_xlim(-0.05, 0.35); ax.set_ylim(-0.15, 0.25); ax.set_zlim(0, 0.4)
fig.tight_layout()
fig.savefig(FIG / "fig2_ik_secuencia_A.png", dpi=150)
print("figura guardada")
