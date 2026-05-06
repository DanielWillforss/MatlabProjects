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

%% a) shortest path

G_l = digraph(edges(:,1), edges(:, 2), l);
P = shortestpath(G_l, 1, 17);
%% b) maxflow

G_C = digraph(edges(:,1), edges(:, 2), C);
mf = maxflow(G_C, 1, 17);

%% c) node net flow

b = B*flow;

%% prep for cvx

M = 28;
nu = zeros(17, 1);
nu(1) = b(1);
nu(17) = -b(1);

cvx_quiet true

%% d)

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
    minimize( sum( l .* C .* (log(C) - log(C - f2))))
    subject to
        B * f2 == nu
        0 <= f2 <= C
cvx_end

%% f)

w1 = f1.*l./(C.*(1-f1./C).^2);

cvx_begin
    variable f3(M)
    minimize( sum( l .* C .* (log(C) - log(C - f3)) + f3.*w1) )
    subject to
        B * f3 == nu
        0 <= f3 <= C
cvx_end

%%

cvx_begin
    variable f4(M)
    minimize sum(l.*quad_over_lin(f4,(C-f4), 0))
    subject to
        B * f4 == nu
        0 <= f4 <= C
cvx_end

%%

w2 = f4.*l./(C.*(1-f4./C).^2) - l;

cvx_begin
    variable f5(M)
    minimize( sum( l .* C .* (log(C) - log(C - f5)) + f5.*w2) )
    subject to
        B * f5 == nu
        0 <= f5 <= C
cvx_end
