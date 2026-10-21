%% MIA-12 Robótica y Sistemas Autónomos · Grupo 5 · Avance del trabajo final
% Comprobación del modelo del robot antropomórfico RRR de 3 GDL de
% Urrea, Kern y Torres (2024), Applied Sciences 14, 11010 (SCORBOT-ER V plus, 3 ejes).
% Son los mismos cinco pasos que los scripts de Python de la carpeta src/.
% Usa solo MATLAB base y Symbolic Math Toolbox.
clear; clc; close all;

L1 = 0.32; L2 = 0.22; L3 = 0.22;                 % Tabla 3 del paper [m]
LIM = [-155 155; -35 130; -130 130];             % Tabla 1 [grados]; supuesto: cero DH = cero del robot
seqA = [0.25 0.01 0.11; 0.25 0.01 0.11; 0.25 0.11 0.11;
        0.15 0.11 0.11; 0.15 0.01 0.11; 0.25 0.01 0.11];   % Tabla 10 (secuencia A) [m]
unicos = seqA([1 3 4 5], :);
nombres = {'P1 = P2 = P6', 'P3', 'P4', 'P5'};

%% Paso 1. Cinemática directa: producto DH frente a las ecuaciones (9), (12), (13) y (14)
fprintf('\n===== PASO 1. CINEMÁTICA DIRECTA =====\n');
fprintf('Qué hacemos: multiplicamos las matrices DH de la Tabla 2 para obtener la posición de la punta del robot\n');
fprintf('(x, y, z) en función de los tres ángulos, y la comparamos con las fórmulas impresas en el paper.\n');
fprintf('Para qué: si una fórmula del paper estuviera mal, todo lo que se construye encima (cinemática inversa,\n');
fprintf('jacobiano, dinámica) saldría mal.\n');
fprintf('\n');
syms t1 t2 t3 l1 l2 l3 real
A1 = dh(t1, l1, 0, sym(pi)/2);
A2 = dh(t2, 0, l2, 0);
A3 = dh(t3, 0, l3, 0);
T01 = A1; T02 = simplify(A1*A2); T03 = simplify(T02*A3);
P = T03(1:3, 4);
fprintf('Px = %s\nPy = %s\nPz = %s\n', char(P(1)), char(P(2)), char(P(3)));

C1 = cos(t1); S1 = sin(t1); C2 = cos(t2); S2 = sin(t2); C3 = cos(t3); S3 = sin(t3);
A2p = [C2 -S2 0 l2*C2; S2 C2 0 l2*C2; 0 0 1 0; 0 0 0 1];   % ec. (9) tal como está impresa
Pxp = l3*C1*C2*C3 - l3*C1*S2*S3 + l2*C1*C2;               % ec. (12)
Pyp = l3*S3*C2*C3 - l3*S1*S2*S3 + l3*S1*C2;               % ec. (13)
Pzp = l3*S2*C3 + l3*C2*S3 + l2*S2 + l1;                   % ec. (14)
fprintf('Ec. (9), fila 2, columna 4: paper - DH = %s -> errata, debe ser l2*sin(theta2)\n', char(simplify(A2p(2,4) - A2(2,4))));
fprintf('Ec. (12) coincide con DH: %s\n', sino(Pxp - P(1)));
fprintf('Ec. (13) coincide con DH: %s\n', sino(Pyp - P(2)));
fprintf('Ec. (14) coincide con DH: %s\n', sino(Pzp - P(3)));
v = [l1 l2 l3 t1 t2 t3]; val = [L1 L2 L3 30*pi/180 20*pi/180 -60*pi/180];
fprintf('Ejemplo q = (30, 20, -60) grados: P_DH = (%.4f, %.4f, %.4f) m; Py con la ec. (13) = %.4f m\n', ...
        double(subs(P, v, val)), double(subs(Pyp, v, val)));

th = linspace(-155, 155, 300)*pi/180;
fdh = matlabFunction(subs(P(2), [l1 l2 l3 t2 t3], [L1 L2 L3 20*pi/180 -60*pi/180]), 'Vars', t1);
fpa = matlabFunction(subs(Pyp, [l1 l2 l3 t2 t3], [L1 L2 L3 20*pi/180 -60*pi/180]), 'Vars', t1);
figure('Name', 'Figura 1');
plot(th*180/pi, fdh(th), 'LineWidth', 2.5); hold on
plot(th*180/pi, fpa(th), '--', 'LineWidth', 2); grid on
xlabel('\theta_1 [°]  (\theta_2 = 20°, \theta_3 = -60°)'); ylabel('P_y [m]');
legend('P_y por producto DH (Tabla 2)', 'P_y con la ec. (13) impresa', 'Location', 'best');
title('Errata de la ec. (13)');

%% Paso 2. Cinemática inversa geométrica sobre los puntos de las Tablas 10 y 11
fprintf('\n===== PASO 2. CINEMÁTICA INVERSA (TABLAS 10 Y 11) =====\n');
fprintf('Qué hacemos: el problema al revés. Para cada punto (x, y, z) de la trayectoria del paper calculamos los\n');
fprintf('ángulos de la base, el hombro y el codo. Hay dos soluciones: codo arriba y codo abajo.\n');
fprintf('Para qué: el robot se mueve girando sus articulaciones. Luego volvemos con la cinemática directa para\n');
fprintf('comprobar que la punta llega al punto (error casi cero) y revisamos los límites de la Tabla 1.\n');
fprintf('\n');
fprintf('%-13s %-18s %9s %9s %9s %13s %8s %10s\n', 'Punto', 'x, y, z [m]', 'th1 [°]', 'th2 [°]', 'th3 [°]', 'error FK [m]', 'límites', 'det J');
etiquetas = {'codo arriba', 'codo abajo'};
for k = 1:size(unicos, 1)
    p = unicos(k, :);
    Q = ik(p, L1, L2, L3);                      % fila 1: codo arriba; fila 2: codo abajo
    for s = 1:2
        qs = Q(s, :); gr = qs*180/pi;
        err = norm(fk(qs, L1, L2, L3) - p);
        ok = all(gr >= LIM(:,1)' & gr <= LIM(:,2)');
        fprintf('%-13s %-18s %9.3f %9.3f %9.3f %13.1e %8s %10.6f  (%s)\n', nombres{k}, ...
                sprintf('%.2f, %.2f, %.2f', p), gr, err, dentro(ok), detJ(qs, L2, L3), etiquetas{s});
    end
end
fprintf('Secuencia B (Tabla 11): igual a la A, salvo P5 = (0.25, 0.01, 0.11), que repite la solución de P1.\n');

figure('Name', 'Figura 2'); hold on; grid on
plot3(seqA(:,1), seqA(:,2), seqA(:,3), 'k--', 'LineWidth', 1.5);
col = lines(4);
for k = 1:4
    Q = ik(unicos(k, :), L1, L2, L3);
    B = puntos_brazo(Q(1, :), L1, L2, L3);
    plot3(B(:,1), B(:,2), B(:,3), '-o', 'LineWidth', 2.5, 'Color', col(k,:), 'MarkerFaceColor', col(k,:));
    text(unicos(k,1), unicos(k,2), unicos(k,3) - 0.03, nombres{k});
end
view(35, 25); axis equal
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Cinemática inversa en la secuencia A (codo arriba)');

%% Paso 3. Espacio de trabajo con los límites de la Tabla 1
fprintf('\n===== PASO 3. ESPACIO DE TRABAJO =====\n');
fprintf('Qué hacemos: movemos las tres articulaciones dentro de sus límites y calculamos todos los puntos que\n');
fprintf('alcanza la punta.\n');
fprintf('Para qué: saber hasta dónde llega el robot y comprobar que los puntos del paper están a su alcance.\n');
fprintf('\n');
n = 45;
[T1g, T2g, T3g] = ndgrid(linspace(LIM(1,1), LIM(1,2), n)*pi/180, ...
                         linspace(LIM(2,1), LIM(2,2), n)*pi/180, ...
                         linspace(LIM(3,1), LIM(3,2), n)*pi/180);
R = L2*cos(T2g) + L3*cos(T2g + T3g);
X = cos(T1g).*R; Y = sin(T1g).*R; Z = L1 + L2*sin(T2g) + L3*sin(T2g + T3g);
rhueco = sqrt(L2^2 + L3^2 + 2*L2*L3*cosd(130));
dist = vecnorm([hypot(unicos(:,1), unicos(:,2)), unicos(:,3) - L1], 2, 2);
fprintf('Alcance máximo desde el hombro: %.4f m (l2 + l3 = %.2f m)\n', max(R(:)), L2 + L3);
fprintf('Radio del hueco interior por el límite del codo (±130°): %.3f m\n', rhueco);
fprintf('Altura del extremo: de %.3f m a %.3f m\n', min(Z(:)), max(Z(:)));
fprintf('Distancia de los puntos de la Tabla 10 al hombro [m]: %s\n', mat2str(dist', 4));

figure('Name', 'Figura 3');
subplot(1, 2, 1);
rng(0); idx = randperm(numel(X), 25000);
scatter3(X(idx), Y(idx), Z(idx), 2, Z(idx), 'filled'); hold on
scatter3(unicos(:,1), unicos(:,2), unicos(:,3), 40, 'r', 'filled');
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Espacio de trabajo 3D');
subplot(1, 2, 2);
rng(1);
t2r = (LIM(2,1) + diff(LIM(2,:))*rand(60000, 1))*pi/180;
t3r = (LIM(3,1) + diff(LIM(3,:))*rand(60000, 1))*pi/180;
scatter(L2*cos(t2r) + L3*cos(t2r + t3r), L1 + L2*sin(t2r) + L3*sin(t2r + t3r), 1, [0.18 0.44 0.69], 'filled'); hold on
scatter(hypot(unicos(:,1), unicos(:,2)), unicos(:,3), 40, 'r', 'filled');
plot(0, L1, 'ks', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('r = (x^2 + y^2)^{1/2} [m]'); ylabel('z [m]');
title('Corte en el plano del brazo');

%% Paso 4. Jacobiano analítico y geométrico; singularidades
fprintf('\n===== PASO 4. JACOBIANO Y SINGULARIDADES =====\n');
fprintf('Qué hacemos: calculamos el jacobiano, la matriz que convierte la velocidad de las articulaciones en\n');
fprintf('velocidad de la punta, por los dos métodos de clase, y buscamos dónde su determinante vale cero.\n');
fprintf('Para qué: donde det(J) = 0 (singularidad) el robot pierde un movimiento y los motores tendrían que\n');
fprintf('girar muy rápido; la trayectoria debe pasar lejos de esas zonas.\n');
fprintf('\n');
q = [t1; t2; t3];
Jan = simplify(jacobian(P, q));                          % derivada de la posición
z0 = [0; 0; 1]; p0 = [0; 0; 0];
z1 = T01(1:3, 3); p1 = T01(1:3, 4);
z2 = T02(1:3, 3); p2 = T02(1:3, 4);
Jv = simplify([cross(z0, P - p0), cross(z1, P - p1), cross(z2, P - p2)]);   % z_{i-1} x (p - p_{i-1})
Jw = [z0 z1 z2];
fprintf('Jacobiano analítico = geométrico (parte lineal): %s\n', sino(Jan - Jv));
disp('Jv ='); disp(Jan);
disp('Jw ='); disp(Jw);
detForma = -l2*l3*sin(t3)*(l2*cos(t2) + l3*cos(t2 + t3));
fprintf('det(Jv) = %s\n', char(simplify(det(Jan))));
fprintf('det(Jv) = -l2*l3*sin(theta3)*(l2*cos(theta2) + l3*cos(theta2 + theta3)): %s\n', sino(det(Jan) - detForma));
fprintf('%-13s %12s %20s\n', 'Punto', 'det(Jv)', 'número de condición');
for k = 1:4
    Q = ik(unicos(k, :), L1, L2, L3);
    J = jac(Q(1, :), L2, L3);
    fprintf('%-13s %12.6f %20.2f\n', nombres{k}, det(J), cond(J));
end

[g2, g3] = ndgrid(linspace(-35, 130, 300), linspace(-130, 130, 300));
D = -L2*L3*sind(g3).*(L2*cosd(g2) + L3*cosd(g2 + g3));
figure('Name', 'Figura 4');
contourf(g2, g3, D, 30, 'LineColor', 'none'); hold on
contour(g2, g3, D, [0 0], 'k', 'LineWidth', 2);
colorbar; colormap(parula);
for k = 1:4
    Q = ik(unicos(k, :), L1, L2, L3); gq = Q(1, :)*180/pi;
    plot(gq(2), gq(3), 'o', 'MarkerFaceColor', [1 0.84 0], 'MarkerEdgeColor', 'k');
end
xlabel('\theta_2 [°]'); ylabel('\theta_3 [°]');
title('det(J_v) y singularidades (línea negra: det = 0)');

%% Paso 5. Revisión del modelo dinámico de Euler-Lagrange (ecs. 21 a 37 del paper)
fprintf('\n===== PASO 5. REVISIÓN DE LA DINÁMICA =====\n');
fprintf('Qué hacemos: armamos las ecuaciones de movimiento de Euler-Lagrange (inercia M, Coriolis C y gravedad G)\n');
fprintf('y las comparamos término a término con las del paper.\n');
fprintf('Para qué: con este modelo se calcula el torque de cada motor; es la base del control y del dataset de\n');
fprintf('la red neuronal.\n');
fprintf('\n');
syms lc2 lc3 m2 m3 gz I1zz Ixx2 Iyy2 I2zz Ixx3 Iyy3 I3zz positive
syms d1 d2 d3 real
dq = [d1; d2; d3];
pc2 = T02*[lc2 - l2; 0; 0; 1]; pc2 = pc2(1:3);      % centro de masa sobre el eje x de cada eslabón
pc3 = T03*[lc3 - l3; 0; 0; 1]; pc3 = pc3(1:3);
Jv2 = jacobian(pc2, q); Jv3 = jacobian(pc3, q);
Jw2 = [z0 z1 sym(zeros(3, 1))]; Jw3 = [z0 z1 z2];
R2 = T02(1:3, 1:3); R3 = T03(1:3, 1:3);
I2 = R2*diag([Ixx2 Iyy2 I2zz])*R2.'; I3 = R3*diag([Ixx3 Iyy3 I3zz])*R3.';
M = sym(zeros(3)); M(1, 1) = I1zz;
M = M + m2*(Jv2.'*Jv2) + Jw2.'*I2*Jw2 + m3*(Jv3.'*Jv3) + Jw3.'*I3*Jw3;
M = simplify(M);
U = gz*(m2*pc2(3) + m3*pc3(3));
G = simplify(jacobian(U, q).');
Cv = sym(zeros(3, 1));
for k = 1:3
    for i = 1:3
        for j = 1:3
            c = (diff(M(k,j), q(i)) + diff(M(k,i), q(j)) - diff(M(i,j), q(k)))/2;   % Christoffel
            Cv(k) = Cv(k) + c*dq(i)*dq(j);
        end
    end
end
Cv = simplify(Cv);

% Modelo impreso en el paper
C23 = cos(t2 + t3); S23 = sin(t2 + t3); r3 = lc3*C23 + l2*C2;
Mp = [r3^2*m3 + I1zz + I2zz + I3zz + (lc2^2*cos(2*t2)/2 + lc2^2/2)*m2, I2zz + I3zz, I3zz;
      I2zz + I3zz, (l2^2 + 2*C3*l2*lc3 + lc3^2)*m3 + m2*lc2^2 + I2zz + I3zz, (lc3^2 + l2*C3*lc3)*m3 + I3zz;
      I3zz, (lc3^2 + l2*C3*lc3)*m3 + I3zz, lc3^2*m3 + I3zz];
Cp = [-d1*d2*(m2*sin(2*t2)*lc2^2 + 2*m3*r3*(lc3*S23 + l2*S2)) - 2*d1*d3*lc3*m3*S23*r3;
      (m2*sin(2*t2)*lc2^2/2 + m3*r3*(lc3*S23 + l2*S2))*d1^2 - l2*lc3*m3*S3*d3^2 - 2*d2*l2*lc3*m3*S3*d3;
      lc3*m3*S23*r3*d1^2 + l2*lc3*m3*S3*d2^2];
Gp = [0; gz*(m3*r3 + lc2*m2*C2); gz*lc3*m3*C23];

% El M11 del paper equivale a tomar Ixx = Iyy = Izz en los eslabones 2 y 3
iso = [Ixx2 Iyy2 Ixx3 Iyy3]; isov = [I2zz I2zz I3zz I3zz];
dM = simplify(subs(M, iso, isov) - Mp);
dC = simplify(subs(Cv, iso, isov) - Cp);
dG = simplify(G - Gp);
disp('M (nuestra derivación, tensor de inercia diagonal general):'); disp(M);
disp('M nuestra (Ixx = Iyy = Izz) - M paper ='); disp(dM);
disp('C nuestra (Ixx = Iyy = Izz) - C paper ='); disp(dC.');
disp('G nuestra - G paper ='); disp(dG.');

% Comprobación numérica en valores al azar (independiente de simplify)
vars = [t1 t2 t3 d1 d2 d3 l2 l3 lc2 lc3 m2 m3 gz I1zz I2zz I3zz];
rng(7); rv = 0.1 + rand(1, numel(vars));
dMn = double(subs(dM, vars, rv)); dCn = double(subs(dC, vars, rv)); dGn = double(subs(dG, vars, rv));
mascara = true(3); mascara(1,2) = false; mascara(1,3) = false; mascara(2,1) = false; mascara(3,1) = false;
fprintf('Comprobación numérica: max|dM| fuera de M12 y M13 = %.1e; max|dC| = %.1e; max|dG| = %.1e\n', ...
        max(abs(dMn(mascara))), max(abs(dCn)), max(abs(dGn)));

% Valores en P1 (codo arriba) con la Tabla 3
Q = ik(seqA(1, :), L1, L2, L3); qP1 = Q(1, :);
nv = [l2 l3 lc2 lc3 m2 m3 I1zz I2zz I3zz gz t1 t2 t3];
nval = [L2 L3 0.052 0.1376 2.2483 1.957 0.04624 0.02595 0.03616 9.81 qP1];
Mn = double(subs(subs(M, iso, isov), nv, nval));
Mpn = double(subs(Mp, nv, nval));
disp('M nuestra en P1 [kg·m^2]:'); disp(round(Mn, 5));
disp('M paper en P1 [kg·m^2]:'); disp(round(Mpn, 5));
fprintf('Torque de inercia para qdd = (0, 1, 0) rad/s^2: nuestro (%.5f, %.5f, %.5f) N·m; paper (%.5f, %.5f, %.5f) N·m\n', ...
        Mn*[0; 1; 0], Mpn*[0; 1; 0]);
fprintf('\nListo: pasos 1 a 5 ejecutados.\n');

%% Funciones locales
function A = dh(theta, d, a, alpha)
% Matriz de transformación homogénea de Denavit-Hartenberg
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
     0           sin(alpha)             cos(alpha)            d;
     0           0                      0                     1];
end

function p = fk(q, L1, L2, L3)
A = dh(q(1), L1, 0, pi/2)*dh(q(2), 0, L2, 0)*dh(q(3), 0, L3, 0);
p = A(1:3, 4)';
end

function Q = ik(p, L1, L2, L3)
% Cinemática inversa geométrica; fila 1: codo arriba, fila 2: codo abajo
x = p(1); y = p(2); z = p(3);
t1 = atan2(y, x); r = hypot(x, y); s = z - L1;
c3 = (r^2 + s^2 - L2^2 - L3^2)/(2*L2*L3);
if abs(c3) > 1, error('Punto fuera del espacio de trabajo'); end
s3v = [-sqrt(1 - c3^2), sqrt(1 - c3^2)];
Q = zeros(2, 3);
for k = 1:2
    t3 = atan2(s3v(k), c3);
    t2 = atan2(s, r) - atan2(L3*s3v(k), L2 + L3*c3);
    Q(k, :) = [t1 t2 t3];
end
end

function J = jac(q, L2, L3)
% Jacobiano lineal 3x3 analítico
c1 = cos(q(1)); s1 = sin(q(1));
r = L2*cos(q(2)) + L3*cos(q(2) + q(3));
dr2 = -L2*sin(q(2)) - L3*sin(q(2) + q(3));
dr3 = -L3*sin(q(2) + q(3));
J = [-s1*r, c1*dr2, c1*dr3;
      c1*r, s1*dr2, s1*dr3;
      0, L2*cos(q(2)) + L3*cos(q(2) + q(3)), L3*cos(q(2) + q(3))];
end

function d = detJ(q, L2, L3)
d = -L2*L3*sin(q(3))*(L2*cos(q(2)) + L3*cos(q(2) + q(3)));
end

function B = puntos_brazo(q, L1, L2, L3)
A1 = dh(q(1), L1, 0, pi/2); T2 = A1*dh(q(2), 0, L2, 0); T3 = T2*dh(q(3), 0, L3, 0);
B = [0 0 0; A1(1:3, 4)'; T2(1:3, 4)'; T3(1:3, 4)'];
end

function s = sino(E)
if isequal(simplify(E), sym(zeros(size(E)))), s = 'sí'; else, s = 'no'; end
end

function s = dentro(b)
if b, s = 'dentro'; else, s = 'fuera'; end
end
