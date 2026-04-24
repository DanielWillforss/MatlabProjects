clc
close all
clear

%% Handin 2

C = load('capacities.mat', '-ascii');
flow = load('flow.mat', '-ascii');
B = load('traffic.mat', '-ascii');
l = load("traveltime.mat", '-ascii'); % traveltime in hours
d = @(f) (l./(1-f./C)).*(f < C) + realmax.*(f > C);

[from_nodes, edge_idx_from] = find(B == 1);
[to_nodes,   edge_idx_to]   = find(B == -1);

edges = zeros(size(B,2), 2);
edges(edge_idx_from, 1) = from_nodes;
edges(edge_idx_to,   2) = to_nodes;

coord = [0 3; 1 3; 2 3; 3 3; 4 3; 0 2; 1 2; 2 2; 3 2; 1 1; 2 1; 3 1; 4 1; 5 1; 2 0; 3 0; 5 0];

%% a) shortest path

G = digraph(edges(:,1), edges(:, 2), l);
P = shortestpath(G, 1, 17);
%% b) maxflow

G = digraph(edges(:,1), edges(:, 2), C);
mf = maxflow(G, 1, 17);

%% c) node net flow

b = zeros(17,1);

% subtract outgoing flow
for i = 1:28
    u = edges(i,1);
    v = edges(i,2);
    
    b(u) = b(u) - flow(i);   % flow leaves u
    b(v) = b(v) + flow(i);   % flow enters v
end

%% d)
%net_in = b(1);

M = 28;
nu = zeros(17, 1);
nu(1) = -b(1);
nu(17) = b(1);

cvx_quiet true

cvx_begin
    variable f1(M)
    minimize sum(l.*C.*inv_pos(1-f1./C) - l.*C)
    subject to
        B * f1 == nu
        0 <= f1 <= C
cvx_end

%% e)

cvx_begin
    variable f2(M)
    minimize( sum( -l .* C .* log(1 - f2./C) ) )
    subject to
        B * f2 == nu
        0 <= f2 <= C
cvx_end

%% f)

w = f1.*l./(C.*(1-f1./C).^2);

cvx_begin
    variable f3(M)
    minimize( sum( -l .* C .* log(1 - f3./C) + f3.*w) )
    subject to
        B * f3 == nu
        0 <= f3 <= C
cvx_end

disp(f3);

%%

cvx_begin
    variable f4(M)
    minimize sum(l.*C.*inv_pos(1-f4./C) - l.*C)
    subject to
        B * f4 == nu
        0 <= f4 <= C
cvx_end

%% Plot capacity

figure; hold on; axis equal;

% Plot nodes as large circles
scatter(coord(:,1), coord(:,2), 800, 'o', ...
    'MarkerFaceColor', [0.8 0 0], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

% Add labels centered inside the circles
for i = 1:size(coord,1)
    text(coord(i,1), coord(i,2), num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

C = C(:); % ensure column vector

Cmin = min(C);
Cmax = max(C);

C_norm = (C - Cmin) / (Cmax - Cmin + eps); % scale to [0,1]

cmap = turbo(256);   % or jet(256), turbo(256), parula(256) etc.

colors = cmap(round(C_norm * 255) + 1, :);

for e = 1:size(edges,1)
    i = edges(e,1); % from node
    j = edges(e,2); % to node

    dx = coord(j,1) - coord(i,1);
    dy = coord(j,2) - coord(i,2);

    quiver(coord(i,1), coord(i,2), dx, dy, 0, ...
        'Color', colors(e,:), ...
        'LineWidth', 1.2, ...
        'MaxHeadSize', 0.2);
end

scatter(coord(edges(:,1),1), coord(edges(:,1),2), ...
    1, f, 'filled', 'MarkerFaceAlpha', 0); % invisible but carries data

colormap(turbo);      % or parula, cividis, etc.
clim([min(C) max(C)]);
colorbar;

% Optional: nicer limits
xlim([min(coord(:,1))-0.5, max(coord(:,1))+0.5]);
ylim([min(coord(:,2))-0.5, max(coord(:,2))+0.5]);

title("Capacity")

axis off

%% Plot graph

figure; hold on; axis equal;

% Plot nodes as large circles
scatter(coord(:,1), coord(:,2), 800, 'o', ...
    'MarkerFaceColor', [0.8 0 0], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

% Add labels centered inside the circles
for i = 1:size(coord,1)
    text(coord(i,1), coord(i,2), num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

f1 = f1(:); % ensure column vector

fmin = min(f1);
fmax = max(f1);

f1_norm = (f1 - fmin) / (f1max - f1min + eps); % scale to [0,1]

cmap = turbo(256);   % or jet(256), turbo(256), parula(256) etc.

colors = cmap(round(f1_norm * 255) + 1, :);

for e = 1:size(edges,1)
    i = edges(e,1); % from node
    j = edges(e,2); % to node

    dx = coord(j,1) - coord(i,1);
    dy = coord(j,2) - coord(i,2);

    quiver(coord(i,1), coord(i,2), dx, dy, 0, ...
        'Color', colors(e,:), ...
        'LineWidth', 1.2, ...
        'MaxHeadSize', 0.2);
end

scatter(coord(edges(:,1),1), coord(edges(:,1),2), ...
    1, f, 'filled', 'MarkerFaceAlpha', 0); % invisible but carries data

colormap(turbo);      % or parula, cividis, etc.
clim([min(f1) max(f1)]);
colorbar;

% Optional: nicer limits
xlim([min(coord(:,1))-0.5, max(coord(:,1))+0.5]);
ylim([min(coord(:,2))-0.5, max(coord(:,2))+0.5]);

title("min flow")

axis off