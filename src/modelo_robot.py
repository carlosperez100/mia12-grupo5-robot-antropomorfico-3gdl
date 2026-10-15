"""
Modelo cinemático del robot antropomórfico RRR de 3 GDL estudiado en
Urrea, Kern y Torres (2024), Appl. Sci. 14, 11010 (SCORBOT-ER V plus, 3 ejes).

Datos tomados del paper:
  - Tabla 2 (DH): (θ1, l1, 0, 90°), (θ2, 0, l2, 0), (θ3, 0, l3, 0)
  - Tabla 3: l1 = 0.32 m, l2 = 0.22 m, l3 = 0.22 m
  - Tabla 1: rango de base 310°, hombro +130°/−35°, codo ±130°
"""
import numpy as np

L1, L2, L3 = 0.32, 0.22, 0.22

# Límites articulares de la Tabla 1 (grados). Se asume que el cero de la
# tabla DH coincide con el cero del robot (supuesto del grupo).
LIMITES = {
    "theta1": (-155.0, 155.0),
    "theta2": (-35.0, 130.0),
    "theta3": (-130.0, 130.0),
}


def dh(theta, d, a, alpha):
    """Matriz de transformación homogénea de Denavit-Hartenberg."""
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0.0,      sa,       ca,      d],
        [0.0,     0.0,      0.0,    1.0],
    ])


def matrices_dh(q):
    """Devuelve 0A1, 1A2, 2A3 para q = (θ1, θ2, θ3) en radianes."""
    t1, t2, t3 = q
    return (dh(t1, L1, 0.0, np.pi / 2),
            dh(t2, 0.0, L2, 0.0),
            dh(t3, 0.0, L3, 0.0))


def cinematica_directa(q):
    """Posición del extremo por producto de matrices DH."""
    A1, A2, A3 = matrices_dh(q)
    return (A1 @ A2 @ A3)[:3, 3]


def cinematica_directa_cerrada(q):
    """Forma cerrada equivalente: x = C1(l2C2 + l3C23), y = S1(...), z = l1 + l2S2 + l3S23."""
    t1, t2, t3 = q
    r = L2 * np.cos(t2) + L3 * np.cos(t2 + t3)
    return np.array([np.cos(t1) * r, np.sin(t1) * r,
                     L1 + L2 * np.sin(t2) + L3 * np.sin(t2 + t3)])


def cinematica_inversa(p):
    """
    Cinemática inversa geométrica (método visto en clase).
    θ1 = atan2(y, x); en el plano del brazo: r = √(x²+y²), s = z − l1,
    C3 por ley de cosenos, S3 = ±√(1 − C3²), θ2 = atan2(s, r) − atan2(l3 S3, l2 + l3 C3).
    Devuelve las dos soluciones: (codo arriba, codo abajo).
    """
    x, y, z = p
    t1 = np.arctan2(y, x)
    r = np.hypot(x, y)
    s = z - L1
    c3 = (r**2 + s**2 - L2**2 - L3**2) / (2 * L2 * L3)
    if abs(c3) > 1:
        raise ValueError("Punto fuera del espacio de trabajo")
    soluciones = []
    for s3 in (-np.sqrt(1 - c3**2), np.sqrt(1 - c3**2)):
        t3 = np.arctan2(s3, c3)
        t2 = np.arctan2(s, r) - np.arctan2(L3 * s3, L2 + L3 * c3)
        soluciones.append(np.array([t1, t2, t3]))
    return soluciones


def jacobiano(q):
    """Jacobiano lineal 3x3 analítico (derivada de la posición respecto a q)."""
    t1, t2, t3 = q
    c1, s1 = np.cos(t1), np.sin(t1)
    r = L2 * np.cos(t2) + L3 * np.cos(t2 + t3)
    dr2 = -L2 * np.sin(t2) - L3 * np.sin(t2 + t3)
    dr3 = -L3 * np.sin(t2 + t3)
    return np.array([
        [-s1 * r, c1 * dr2, c1 * dr3],
        [ c1 * r, s1 * dr2, s1 * dr3],
        [0.0, L2 * np.cos(t2) + L3 * np.cos(t2 + t3), L3 * np.cos(t2 + t3)],
    ])


def det_jacobiano(q):
    """det(J) = −l2·l3·sin(θ3)·(l2·cos θ2 + l3·cos(θ2+θ3)) (obtenido en el paso 4)."""
    t1, t2, t3 = q
    return -L2 * L3 * np.sin(t3) * (L2 * np.cos(t2) + L3 * np.cos(t2 + t3))


def dentro_de_limites(q):
    g = np.degrees(q)
    return all(lo <= v <= hi for v, (lo, hi) in zip(g, LIMITES.values()))
