%% Hello world

clear

%static variables
Q = [0, 2/5, 1/5, 0, 0;
         0, 0, 3/4, 1/4, 0;
         1/2, 0, 0, 1/2, 0;
         0, 0, 1/3, 0, 2/3;
         0, 1/3, 0, 1/3, 0];

w = Q*ones(5,1);
P = Q./w;
cumulative = cumsum(P')';

% config
init_node=2;
end_node= 2;
N = 1e5;
simulation_limit = 1000;

% setup
T = zeros(N, 1);
for n = 1:N

    % reset simulation
    current_node = init_node;
    t = 0; 
    
    for i = 1:simulation_limit
        t_next = -log(rand())/w(current_node);
        next_state = find(cumulative(current_node,:) >= rand(), 1);
        
        t = t + t_next;
        current_node = next_state;
    
        if current_node == end_node
            break;
        end
    end
    T(n) = t;

end

avg_time = mean(T);
std_dev = std(T);
stderr = std_dev / sqrt(N);

disp("avg: " + avg_time + " stderr: " + stderr)

%% Calculate mean

L = diag(w) - Q;
dist = null(L');
dist_norm = dist / sum(dist);
avg_return = 1./(dist_norm.*w);

target = 5;
idx = setdiff(1:size(P,1), target);
P_sub = P(idx, idx);
w_subinv = 1./w(idx);
t = (eye(size(P_sub)) - P_sub) \ w_subinv;


%%

n = 10;
W = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);

noise = @(t) t/100;
c = @(x, y) double(x == y);
U = @(X) sum(W.*c(X, X'), 'all');
p_help = @(X, i, t, s) exp(-noise(t)*W(i, :)*c(s, X));
P_0 = @(X, i, t) p_help(X, i, t, 0)/sum(p_help(X, i, t, 0) + p_help(X, i, t, 1));

X = zeros(n, 1);
t = 0;
u = zeros(1000, 1);

figure

for i = 1:1000

    id = randi([1, n]);
    P = P_0(X, id, t);
    if rand() < P
        X(id) = 0;
    else
        X(id) = 1;
    end

    t = t + 1;
    scatter(1:n, X);
    ylim([-0.1 1.1]);
    drawnow

    u(i) = U(X);
    if 0.0001 > u(i)
        u = u(1:i);
        break
    end

end

figure
plot(u);

%%

n = 100;
W = load("data/wifi.mat", '-ascii');

noise = @(t) 200;
c = @(x, y) 2*(x == y) + 1*(abs(x - y) == 1);
U = @(X) sum(W.*c(X, X'), 'all');
p_help = @(X, i, t, s) exp(-noise(t)*W(i, :)*c(s, X));

X = ones(n, 1);
t = 0;
u = zeros(1000, 1);

figure

for i = 1:1000

    id = randi([1, n]);
    
    P = zeros(8, 1);
    for j=1:8
        P(j) = p_help(X, id, t, j);
    end
    P = P/sum(P);

    X(id) = find(cumsum(P) >= rand(), 1);

    t = t + 1;
    %scatter(1:n, X);
    %ylim([0.5 8.5]);
    %drawnow

    u(i) = U(X);
    if 0.0001 > u(i)
        u = u(1:i);
        break
    end

end

figure
plot(u);

%%

pos = load("data/coord.mat", '-ascii');


n = length(X);
% Create graph object
G = graph(W);

% Define 8 colors
cmap = hsv(8);   % or hsv(8), jet(8), etc.

figure;
hold on;
axis equal;

% Plot edges first
[row,col] = find(triu(W));   % only upper triangle to avoid duplicates

for k = 1:length(row)
    i = row(k);
    j = col(k);

    plot([pos(i,1) pos(j,1)], ...
         [pos(i,2) pos(j,2)], ...
         'k-', 'LineWidth', 1);
end

% Plot nodes by color group
for c = 1:8
    idx = (X == c);

    scatter(pos(idx,1), ...
            pos(idx,2), ...
            80, ...
            cmap(c,:), ...
            'filled');
end

hold off;