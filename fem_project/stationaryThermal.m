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
ndof = length(p);
coord = p';

edof = zeros(nelm, 4);
edof(:,1) = 1:nelm;
edof(:,2:4) = t(1:3, :)';

dof = (1:ndof)';

nen = 3;
[ex, ey] = coordxtr(edof, coord, dof, 3);

% Set material constants

k = 1;
D = [k 0; 0 k];
qout = 2000;
alphac = 50;
Tinf = 293;
ep = 1;
Q = 4*10^5;

% Set up K, f, and bc
K = zeros(ndof);
f = zeros(ndof, 1);

for enr = 1:length(e)
    boundary = e(5, enr);
    node1 = e(1, enr);
    node2 = e(2, enr);
    Le = norm(coord(node1, :) - coord(node2, :));
    Ke = zeros(2, 2);
    fe = zeros(2, 1);

    if boundary == 3
        fe = -qout*Le/2 * ones(2, 1);
    end
    if ismember(boundary, [1 5:20])
        Ke = alphac*Le/6*[2 1; 1 2];
        fe = alphac*Tinf*Le/2*ones(2,1);
    end
    [K, f] = assem([0 node1 node2], K, Ke, f, fe);
end

%% Solver

for elnr = 1:nelm
    [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), ep, [k 0; 0 k], Q);
    [K, f] = assem(edof(elnr,:), K, Ke, f, fe);
end

a = solve(K, f);

%% Post-processor

ed = extract(edof, a);

%plot
figure
patch(ex', ey', ed', ed', 'EdgeColor', 'none');

axis equal
colorbar
colormap(jet)

xlim([min(ex(:)) max(ex(:))])
ylim([min(ey(:)) max(ey(:))])

xlabel('x')
ylabel('y')
title('Temperature Field')

set(gca, 'FontSize', 12)
