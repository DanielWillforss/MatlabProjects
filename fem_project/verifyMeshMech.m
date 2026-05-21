% Mechanical
clear
close all
clc

%% Pre-processor

%mesh
p = [0 0.005 0.01 0 0.005 0.01 0 0.005 0.01; 
     0 0 0 0.005 0.005 0.005 0.01 0.01 0.01];
e = [7 8 9 6 3 2 1 4;
     8 9 6 3 2 1 4 7;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0;
     1 2 3 4 5 6 7 8];
t = [1 2 3 6 9 8 7 4;
     2 3 6 9 8 7 4 1;
     5 5 5 5 5 5 5 5];


% Plot of Mesh
figure
patch('Faces', t(1:3,:)', ...
      'Vertices', p', ...
      'FaceColor', [0.8 0.9 1.0], ...
      'EdgeColor', [0.2 0.2 0.2]);
axis equal
title('Mesh')

% Compute mesh variables
nelm = length(t);
ndof = 2*length(p);
coord = p';

edof = zeros(nelm, 7);
conn = t(1:3,:)';

edof(:,1) = 1:nelm;
edof(:,2:7) = [
    2*conn(:,1)-1, 2*conn(:,1), ...
    2*conn(:,2)-1, 2*conn(:,2), ...
    2*conn(:,3)-1, 2*conn(:,3)
];
dof = zeros(ndof/2, 2);
dof(:, 1) = ((1:ndof/2)*2-1)';
dof(:, 2) = ((1:ndof/2)*2)';

nen = 3; %triangle
[ex, ey] = coordxtr(edof, coord, dof, nen);

% Set material constants



thickness = 0.01;
E = 200*1e+09;
v = 0.3;

D = hooke(2, E, v);

% Set up K, f, and bc
K = zeros(ndof);
f = zeros(ndof, 1);
bc = [];

figure
title("normals")
hold on
axis equal


for enr = 1:length(e)
    boundary = e(5, enr);
    node1 = e(1, enr);
    node2 = e(2, enr);
    Le = norm(coord(node1, :) - coord(node2, :));
    
    dx = coord(node2, 1) - coord(node1, 1);
    dy = coord(node2, 2) - coord(node1, 2);
    normal = [-dy dx];
    normal = normal / norm(normal, 2);

    %fe = thickness*Le/2 * [-pnx -pny -pnx -pny]
    pn = 0; %securded edge

    if ismember(boundary, [2])
        pn = 150*1e6;
    end
    if ismember(boundary, [3 4])
        pn = 100*1e6; 
    end
    if ismember(boundary, [7 8])
        bc = [bc;
              2*node1-1  0;
              2*node1    0;
              2*node2-1  0;
              2*node2    0];
    end

    fe = thickness*Le/2*pn * [normal(1); normal(2); normal(1); normal(2)]; 
    
    edof_edge = [
            2*node1-1
            2*node1
            2*node2-1
            2*node2
        ];
    
    f(edof_edge) = f(edof_edge) + fe;
    
    %plot
    midpoint = [coord(node1, 1) + dx*0.5 coord(node1, 2) + dy*0.5];
    %quiver(midpoint(1), midpoint(2), normal(1), normal(2), 0.001, 'MaxHeadSize', 2);
    quiver(midpoint(1), midpoint(2), fe(1), fe(2), 1);
end


hold off

%% Thermal
%load("element_temp.mat");

%T = sum(ed, 2)/3;
%T0 = 293;
%alpha = 8e-6;

eps_th = 0;
%eps_th = alpha.*(T - T0)*[1 1 1 0];



%% Solver

for elnr = 1:nelm
    Ke = plante(ex(elnr,:), ey(elnr,:), [2 thickness], D);
    %ef = plantf(ex(elnr,:), ey(elnr,:), [2 thickness], (D*eps_th(elnr,:)')');
    %[K, f] = assem(edof(elnr,:), K, Ke, f, ef);
    K = assem(edof(elnr,:), K, Ke);

   % [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), ep, [k 0; 0 k], Q);
    %[K, f] = assem(edof(elnr,:), K, Ke, f, fe);
end

[u, Q] = solveq(K, f, bc);

%% Post-processor

ed = extract(edof, u);
[~, eps] = plants(ex, ey, [2 thickness], D, ed);
sigma = (D*(eps-eps_th)')';

xx2 = sigma(:, 1).^2;
yy2 = sigma(:, 2).^2;
zz2 = sigma(:, 3).^2;
xxyy = sigma(:, 1).*sigma(:, 2);
xxzz = sigma(:, 1).*sigma(:, 3);
yyzz = sigma(:, 2).*sigma(:, 3);
xy2 = sigma(:, 4).^2;

sigma_eff = sqrt(xx2 + yy2 + zz2 - xxyy - xxzz - yyzz + 3*xy2);

%%
% Extract nodal displacements
ux = u(1:2:end);
uy = u(2:2:end);

% Displacement magnitude
umag = sqrt(ux.^2 + uy.^2);

% Plot
figure

patch('Faces', t(1:3,:)', ...
      'Vertices', p', ...
      'FaceVertexCData', umag, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'none');

hold on

quiver(p(1,:), p(2,:), ux', uy', 'k')

axis equal
colorbar
colormap(jet)

xlabel('x')
ylabel('y')
title('Displacement magnitude and direction')

% Deformation scale factor
scale = 100;   % adjust as needed

% Original coordinates
x = p(1,:)';
y = p(2,:)';

% Deformed coordinates
xd = x + scale*ux;
yd = y + scale*uy;

figure
hold on

% Undeformed mesh (gray wireframe)
patch('Faces', t(1:3,:)', ...
      'Vertices', [x y], ...
      'FaceColor', 'none', ...
      'EdgeColor', [0.7 0.7 0.7]);

% Deformed mesh
patch('Faces', t(1:3,:)', ...
      'Vertices', [xd yd], ...
      'FaceVertexCData', umag, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'k');

axis equal
colorbar
colormap(jet)

xlabel('x')
ylabel('y')

title(['Deformed shape (scale = ' num2str(scale) ')'])

figure
hold on

% Undeformed mesh (gray wireframe)
patch('Faces', t(1:3,:)', ...
      'Vertices', [x y], ...
      'FaceColor', 'none', ...
      'EdgeColor', [0.7 0.7 0.7]);

% Deformed mesh with von Mises stress
patch('Faces', t(1:3,:)', ...
      'Vertices', [xd yd], ...
      'FaceVertexCData', sigma_eff, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'k');

axis equal
colorbar
colormap(jet)

xlabel('x')
ylabel('y')

title('Von Mises stress')