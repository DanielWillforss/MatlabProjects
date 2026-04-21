clear
close all
clc

%% Pre-processor (old)

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

% Set up K, f, and bc
topBoundary = e(1:2, (e(5,:) == 1));
topNodes = unique(topBoundary);
bottomBoundary = e(1:2, (e(5,:) == 3));
bottomNodes = unique(bottomBoundary);
circleBoundary = e(1:2, ismember(e(5,:), 5:20));
circleNodes = unique(circleBoundary);

bc = [topNodes, 1000*ones(length(topNodes), 1); 
    bottomNodes, zeros(length(bottomNodes), 1);
    circleNodes, 500*ones(length(circleNodes), 1)];

K = zeros(ndof);
F = zeros(ndof, 1);

%% Solver (old)

for elnr = 1:nelm
    Ke = flw2te(ex(elnr,:), ey(elnr,:), 1, [k 0; 0 k]);
    K = assem(edof(elnr,:), K, Ke);
end

a = solve(K, F, bc);

%% Post-processor (old)

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
