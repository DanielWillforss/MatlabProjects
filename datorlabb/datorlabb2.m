clear
clc

%% Pre-processor

%get mesh
geom2;

k = 1;

K = zeros(ndof);
F = zeros(ndof, 1);

%% Solver

for elnr = 1:nelm
    Ke = flw2te(ex(elnr,:), ey(elnr,:), 1, [k 0; 0 k]);
    K = assem(edof(elnr,:), K, Ke);
end

a = solve(K, F, bc);

%% Post-processor

ed = extract(edof, a);

%Restoring full shape from symmetry
ex = [ex; ex; -ex; -ex];
ey = [ey; -ey; -ey; ey];
ed = [ed; ed; ed; ed];

%plot
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
