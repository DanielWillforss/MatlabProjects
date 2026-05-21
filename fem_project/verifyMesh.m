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
ndof = length(p);
coord = p';

edof = zeros(nelm, 4);
edof(:,1) = 1:nelm;
edof(:,2:4) = t(1:3, :)';

dof = (1:ndof)';

nen = 3;
[ex, ey] = coordxtr(edof, coord, dof, 3);

% Set material constants

k = 10; %k = 2?
D = [k 0; 0 k];
%qout = 2000;
alphac = 1000;
Tinf = 20;
ep = 0.01;
Q = 0;

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

    if ismember(boundary, [3 4])
        bc = [bc;
              node1  100;
              node2  100];
    end
    if ismember(boundary, [2])
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

%not solveq?
a = solve(K, f, bc);

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
