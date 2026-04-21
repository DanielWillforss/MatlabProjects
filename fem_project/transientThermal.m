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
alphac = 50;
Tinf = 293;
T0 = 293;
ep = 1;
Qmax = 4*10^5;
qmax = 2000;
t_tot = 600; %minutes
d0 = ones(ndof, 1)*T0;
Q = @(t) Qmax*sin(pi*t/t_tot);
qout = @(t) qmax*sin(pi*t/t_tot);
nsteps = 6;
t = t_tot/nsteps:t_tot/nsteps:t_tot;

p = 1500;
cp = 800;

% Set up K, f, and bc
K = zeros(ndof);
f = zeros(ndof, nsteps + 1);
M = zeros(ndof);

%% Solver

for elnr = 1:nelm
    Ke = flw2te(ex(elnr,:), ey(elnr,:), ep, [k 0; 0 k]); 
    K = assem(edof(elnr,:), K, Ke);
    Me = plantml(ex(elnr,:), ey(elnr,:), p*cp);
    M = assem(edof(elnr, :), M, Me);
end

for enr = 1:length(e)
    boundary = e(5, enr);
    node1 = e(1, enr);
    node2 = e(2, enr);
    Le = norm(coord(node1, :) - coord(node2, :));
    Ke = zeros(2, 2);
    fe = zeros(2, 1);

    if ismember(boundary, [1 5:20])
        Ke = alphac*Le/6*[2 1; 1 2];
        fe = alphac*Tinf*Le/2*ones(2,1);
    end
    [K, f(:, 1)] = assem([0 node1 node2], K, Ke, f(:, 1), fe);
end



for step = 1:nsteps
    for elnr = 1:nelm
        [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), ep, [k 0; 0 k], Q(t(step))); 
        [~, f(:, step+1)] = assem(edof(elnr,:), K, Ke, f(:, step+1), fe);
    end

    for enr = 1:length(e)
        boundary = e(5, enr);
        node1 = e(1, enr);
        node2 = e(2, enr);
        Le = norm(coord(node1, :) - coord(node2, :));
        Ke = zeros(2, 2);
        fe = zeros(2, 1);
    
        if boundary == 3
            fe = -qout(t(step))*Le/2 * ones(2, 1);
        end
        if ismember(boundary, [1 5:20])
            fe = alphac*Tinf*Le/2*ones(2,1);
        end
        [~, f(:, step+1)] = assem([0 node1 node2], K, Ke, f(:, step+1), fe);
    end

end

ip = [t(1) t_tot 1 [nsteps 0 t]]; %timestep parameters

Tsnap=step1(K,M,d0,ip,f,[]); %pbound?
%%a = solve(K, f);

%% Post-processor

for i = 1:6
ed = extract(edof, Tsnap(:, i));

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
end