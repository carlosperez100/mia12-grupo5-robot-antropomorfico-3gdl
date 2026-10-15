"""
Paso 5. Revisión del modelo dinámico de Euler-Lagrange publicado en el paper (ecs. 20-37).

Se arma el lagrangiano con la geometría DH de la Tabla 2:
  - centro de masa del eslabón i sobre su eje x_i, a lc_i de la articulación i (como en la Tabla 3);
  - tensor de inercia de los eslabones 2 y 3 diagonal en su marco DH (Ixx, Iyy, Izz);
  - eslabón 1: solo gira alrededor de z0, inercia I1zz.
M = Σ m_i Jv_ciᵀ Jv_ci + Jw_iᵀ R_i I_i R_iᵀ Jw_i,  G = ∂U/∂q,  C por símbolos de Christoffel.
Luego se compara término a término con lo impreso en el paper.
"""
import json
import numpy as np
import sympy as sp
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "resultados"

t1, t2, t3 = sp.symbols("theta1 theta2 theta3", real=True)
d1, d2, d3 = sp.symbols("dtheta1 dtheta2 dtheta3", real=True)
l1, l2, l3, lc2, lc3, m2, m3, g = sp.symbols("l1 l2 l3 lc2 lc3 m2 m3 g", positive=True)
I1zz, Ixx2, Iyy2, Izz2, Ixx3, Iyy3, Izz3 = sp.symbols("I1zz Ixx2 Iyy2 I2zz Ixx3 Iyy3 I3zz", positive=True)
q, dq = sp.Matrix([t1, t2, t3]), sp.Matrix([d1, d2, d3])


def dh(theta, d, a, alpha):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta) * sp.cos(alpha),  sp.sin(theta) * sp.sin(alpha), a * sp.cos(theta)],
        [sp.sin(theta),  sp.cos(theta) * sp.cos(alpha), -sp.cos(theta) * sp.sin(alpha), a * sp.sin(theta)],
        [0, sp.sin(alpha), sp.cos(alpha), d], [0, 0, 0, 1]])


T01 = dh(t1, l1, 0, sp.pi / 2)
T02 = sp.simplify(T01 * dh(t2, 0, l2, 0))
T03 = sp.simplify(T02 * dh(t3, 0, l3, 0))

# Centros de masa (en el marco DH del eslabón, el origen está al final del eslabón)
pc2 = (T02 * sp.Matrix([lc2 - l2, 0, 0, 1]))[:3, 0]
pc3 = (T03 * sp.Matrix([lc3 - l3, 0, 0, 1]))[:3, 0]
Jv2, Jv3 = pc2.jacobian(q), pc3.jacobian(q)

z0, z1, z2 = sp.Matrix([0, 0, 1]), T01[:3, 2], T02[:3, 2]
Jw2 = sp.Matrix.hstack(z0, z1, sp.zeros(3, 1))
Jw3 = sp.Matrix.hstack(z0, z1, z2)
I2 = T02[:3, :3] * sp.diag(Ixx2, Iyy2, Izz2) * T02[:3, :3].T
I3 = T03[:3, :3] * sp.diag(Ixx3, Iyy3, Izz3) * T03[:3, :3].T

M = sp.zeros(3, 3)
M[0, 0] += I1zz
M += m2 * Jv2.T * Jv2 + Jw2.T * I2 * Jw2
M += m3 * Jv3.T * Jv3 + Jw3.T * I3 * Jw3
M = sp.simplify(M)

U = g * (m2 * pc2[2] + m3 * pc3[2])
G = sp.simplify(sp.Matrix([sp.diff(U, qi) for qi in q]))

# Coriolis y centrífuga: C_k = Σ_ij c_ijk dq_i dq_j (Christoffel)
C = sp.zeros(3, 1)
for k in range(3):
    for i in range(3):
        for j in range(3):
            c = (sp.diff(M[k, j], q[i]) + sp.diff(M[k, i], q[j]) - sp.diff(M[i, j], q[k])) / 2
            C[k] += c * dq[i] * dq[j]
C = sp.simplify(C)

# --- Modelo impreso en el paper (ecs. 21-29, 31-33, 35-37)
C2_, S2_, C3_, S3_ = sp.cos(t2), sp.sin(t2), sp.cos(t3), sp.sin(t3)
C23, S23 = sp.cos(t2 + t3), sp.sin(t2 + t3)
M_p = sp.Matrix([
    [(lc3 * C23 + l2 * C2_)**2 * m3 + I1zz + Izz2 + Izz3 + (lc2**2 * sp.cos(2 * t2) / 2 + lc2**2 / 2) * m2,
     Izz2 + Izz3, Izz3],
    [Izz2 + Izz3, (l2**2 + 2 * C3_ * l2 * lc3 + lc3**2) * m3 + m2 * lc2**2 + Izz2 + Izz3,
     (lc3**2 + l2 * C3_ * lc3) * m3 + Izz3],
    [Izz3, (lc3**2 + l2 * C3_ * lc3) * m3 + Izz3, lc3**2 * m3 + Izz3]])
r3 = lc3 * C23 + l2 * C2_
C_p = sp.Matrix([
    -d1 * d2 * (m2 * sp.sin(2 * t2) * lc2**2 + 2 * m3 * r3 * (lc3 * S23 + l2 * S2_)) - 2 * d1 * d3 * lc3 * m3 * S23 * r3,
    (m2 * sp.sin(2 * t2) * lc2**2 / 2 + m3 * r3 * (lc3 * S23 + l2 * S2_)) * d1**2
    - l2 * lc3 * m3 * S3_ * d3**2 - 2 * d2 * l2 * lc3 * m3 * S3_ * d3,
    lc3 * m3 * S23 * r3 * d1**2 + l2 * lc3 * m3 * S3_ * d2**2])
G_p = sp.Matrix([0, g * (m3 * r3 + lc2 * m2 * C2_), g * lc3 * m3 * C23])

# Comparación con la simplificación que el paper usa para M11: Ixx = Iyy = Izz en los eslabones 2 y 3
iso = {Ixx2: Izz2, Iyy2: Izz2, Ixx3: Izz3, Iyy3: Izz3}
dM = sp.simplify(M.subs(iso) - M_p)
dC = sp.simplify(C.subs(iso) - C_p)
dG = sp.simplify(G - G_p)

print("M (nuestro, tensor diagonal general):")
for i in range(3):
    print("  ", [sp.simplify(M[i, j]) for j in range(3)])
print("M_nuestro(Ixx=Iyy=Izz) - M_paper =", dM.tolist())
print("C_nuestro(Ixx=Iyy=Izz) - C_paper =", dC.T.tolist())
print("G_nuestro - G_paper =", dG.T.tolist())

# --- Evaluación numérica en P1 (codo arriba) con la Tabla 3
from modelo_robot import cinematica_inversa
qP1 = cinematica_inversa((0.25, 0.01, 0.11))[0]
val = {l1: 0.32, l2: 0.22, l3: 0.22, lc2: 0.052, lc3: 0.1376, m2: 2.2483, m3: 1.957,
       I1zz: 0.04624, Izz2: 0.02595, Izz3: 0.03616, g: 9.81,
       t1: qP1[0], t2: qP1[1], t3: qP1[2]}
Mn = np.array(M.subs(iso).subs(val).evalf(), dtype=float)
Mpn = np.array(M_p.subs(val).evalf(), dtype=float)
print("M nuestro en P1:\n", Mn.round(5))
print("M paper en P1:\n", Mpn.round(5))
# Ejemplo físico: brazo en reposo, se acelera solo el hombro con 1 rad/s²
acc = np.array([0.0, 1.0, 0.0])
tau_n, tau_p = Mn @ acc, Mpn @ acc
print("Torque de inercia para q̈ = (0, 1, 0) rad/s²: nuestro", tau_n.round(5), " paper", tau_p.round(5))

res = {
    "M_nuestro": [[str(sp.simplify(M[i, j])) for j in range(3)] for i in range(3)],
    "G_nuestro": [str(e) for e in G],
    "diferencia_M_con_Ixx_Iyy_Izz_iguales": [[str(e) for e in fila] for fila in dM.tolist()],
    "diferencia_C": [str(e) for e in dC],
    "diferencia_G": [str(e) for e in dG],
    "P1_q_grados": np.degrees(qP1).round(3).tolist(),
    "M_nuestro_P1": Mn.round(5).tolist(),
    "M_paper_P1": Mpn.round(5).tolist(),
    "tau_inercia_qdd_0_1_0_nuestro": tau_n.round(5).tolist(),
    "tau_inercia_qdd_0_1_0_paper": tau_p.round(5).tolist(),
}
(OUT / "paso5_dinamica.json").write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
