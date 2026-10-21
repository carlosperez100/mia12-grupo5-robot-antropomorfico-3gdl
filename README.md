# Comprobación del modelo de un robot antropomórfico RRR de 3 GDL

**Universidad Nacional de Ingeniería · FIIS · Unidad de Posgrado**
MIA-12 Robótica y Sistemas Autónomos · Maestría en Inteligencia Artificial · 2026-2, sección B
Docente: Mg. Wilson Dennis Reyes Vásquez

**Grupo 5:** René Abarca Urbano, Ángel Joel Alvarado Crisanto, Javier Alexander Fernández Purizaca, Carlos Pérez Pérez y Jorge Valega.

Página con el avance y los resultados: **https://carlosperez100.github.io/mia12-grupo5-robot-antropomorfico-3gdl/**

## De qué trata

Avance del trabajo final del curso. Tomamos como paper base:

> Urrea, C.; Kern, J.; Torres, V. (2024). *Design, Simulation, and Comparison of Advanced Control Strategies for a 3-Degree-of-Freedom Robot.* Applied Sciences, 14(23), 11010. https://doi.org/10.3390/app142311010

El paper modela los tres primeros ejes (base, hombro y codo) de un SCORBOT-ER V plus. Con los métodos del curso comprobamos su modelo y completamos lo que no desarrolla.

| Paso | Script | Resultado principal |
|---|---|---|
| 1. Cinemática directa | `src/paso1_verificacion_cinematica_directa.py` | Las ecs. (12) y (14) coinciden con el producto DH. Las ecs. (9) y (13) tienen erratas |
| 2. Cinemática inversa | `src/paso2_cinematica_inversa_trayectorias.py` | Ángulos de los puntos de las Tablas 10 y 11. Solo la solución codo arriba respeta los límites. Error al volver con la directa < 1e-16 m |
| 3. Espacio de trabajo | `src/paso3_espacio_trabajo.py` | Alcance de 0.44 m desde el hombro y hueco interior de 0.186 m. Las trayectorias del paper son alcanzables |
| 4. Jacobiano y singularidades | `src/paso4_jacobiano_singularidades.py` | El jacobiano analítico es igual al geométrico. det(Jv) = −l2·l3·sin θ3·(l2 cos θ2 + l3 cos(θ2+θ3)). Las trayectorias pasan lejos de las singularidades |
| 5. Revisión de la dinámica | `src/paso5_revision_dinamica.py` | M, C y G coinciden con el paper, salvo M12 y M13, que según nuestra derivación son 0 |

Lo que sigue: control por par calculado y PID no lineal, un dataset generado con el modelo, una red neuronal comparada con la solución clásica, la validación con RMS, RSD e IA, y una demo en simulación.

## Cómo ejecutar

```bash
pip install -r requirements.txt
cd src
python paso1_verificacion_cinematica_directa.py
python paso2_cinematica_inversa_trayectorias.py
python paso3_espacio_trabajo.py
python paso4_jacobiano_singularidades.py
python paso5_revision_dinamica.py
```

Las salidas se guardan en `resultados/` (JSON y CSV) y en `assets/figuras/` (PNG). Probado con Python 3.13, numpy 2.3, sympy 1.14 y matplotlib 3.10.

## Estructura

```
index.html          página del avance (GitHub Pages)
src/                modelo del robot y un script por paso
resultados/         salidas numéricas y simbólicas
assets/figuras/     figuras usadas en la página
articulos/          catálogo de los 41 artículos revisados y PDF con licencia abierta
```

En [`articulos/`](articulos/) está el catálogo de los 41 artículos que revisamos para elegir el paper base: 35 de la búsqueda sistemática y 6 de la revisión inicial. Incluye los PDF cuya licencia permite compartirlos, entre ellos el paper base (CC BY 4.0); el resto se cita con su enlace.
