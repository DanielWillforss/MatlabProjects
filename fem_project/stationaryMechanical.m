clear
close all
clc

%% Pre-processor

% Get mesh
load("mesh.mat");

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

tn = 1;
E = 30;
v = 0.3;

D = hooke(2, E, v);

% Set up K, f, and bc
K = zeros(ndof);
f = zeros(ndof, 1);
bc = [];

for enr = 1:length(e)
    boundary = e(5, enr);
    node1 = e(1, enr);
    node2 = e(2, enr);
    Le = norm(coord(node1, :) - coord(node2, :));
    
    
    Ke = zeros(2, 2);
    fe = zeros(2, 1);

    %fe = tn*Le/2 * [-pnx -pny -pnx -pny]

    if boundary == 2
    
        p0 = 1e8;
    
        ymid = (coord(node1,2) + coord(node2,2))/2;
    
        pn = p0*(ymid/0.023 - 1);
    
        tx = 0;
        ty = pn;
    
        fe = tn*Le/2 * [
            tx;
            ty;
            tx;
            ty
        ];
    
        edof_edge = [
            2*node1-1
            2*node1
            2*node2-1
            2*node2
        ];
    
        f(edof_edge) = f(edof_edge) + fe;
    end
    if boundary == 3
        pn = -1e7;     % downward pressure

        tx = 0;
        ty = -pn;
    
        fe = tn*Le/2 * [
            tx;
            ty;
            tx;
            ty
        ];
    
        edof_edge = [
            2*node1-1
            2*node1
            2*node2-1
            2*node2
        ];
    
        f(edof_edge) = f(edof_edge) + fe;
    end
    if boundary == 4
    
        bc = [bc;
              2*node1-1  0;
              2*node1    0;
              2*node2-1  0;
              2*node2    0];
    end
    if ismember(boundary, [5:20])
    
        pn = 1e6;

        tye = coord(node2, 2)-coord(node1, 2);
        txe = coord(node2, 1)-coord(node1, 1);
    
        normal = [tye; -txe];
        normal = normal/norm(normal);
    
        traction = pn*normal;
    
        fe = tn*Le/2 * [
            traction(1);
            traction(2);
            traction(1);
            traction(2)
        ];
    
        edof_edge = [
            2*node1-1
            2*node1
            2*node2-1
            2*node2
        ];
    
        f(edof_edge) = f(edof_edge) + fe;
    end
    
end

%% Solver

for elnr = 1:nelm
    Ke = plante(ex(elnr,:), ey(elnr,:), [2 tn], D);
    K = assem(edof(elnr,:), K, Ke);

   % [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), ep, [k 0; 0 k], Q);
    %[K, f] = assem(edof(elnr,:), K, Ke, f, fe);
end

[a, Q] = solveq(K, f, bc);

%% Post-processor

% Extract nodal displacements
ux = a(1:2:end);
uy = a(2:2:end);

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