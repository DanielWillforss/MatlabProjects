%% Hello world

clear

%static variables
trans = [0, 2/5, 1/5, 0, 0;
         0, 0, 3/4, 1/4, 0;
         1/2, 0, 0, 1/2, 0;
         0, 0, 1/3, 0, 2/3;
         0, 1/3, 0, 1/3, 0];

rate = trans*ones(5,1);
probs = trans./rate;
cumulative = cumsum(probs')';

% config
init_node=2;
end_node= 2;
N = 1e6;
simulation_limit = 1000;

% setup
T = zeros(N, 1);
for n = 1:N

    % reset simulation
    current_node = init_node;
    t = 0; 
    
    for i = 1:simulation_limit
        t_next = -log(rand())/rate(current_node);
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

L = diag(rate) - trans;
dist = null(L');
dist_norm = dist / sum(dist);
avg = 1/(dist_norm(2)*rate(2));