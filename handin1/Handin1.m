%% Functions

function [w_out, w_in] = degree(W)
    [m, n] = size(W); % m must equal n?
    w_out = W * ones(n, 1);
    w_in = W' * ones(m, 1);
end

function out = findLargest(v, n)
    %v must be vector?
    [sortedVals, indices] = sort(v, 'descend');

    topIndices = indices(1:n); % n is larger than length(v)?
    topValues  = sortedVals(1:n);
    out = [topIndices, topValues];
end

function W_comp = largestComponent(W)
    G = digraph(W);
    bins = conncomp(G, 'Type','strong');

    counts = histcounts(bins, 1:(max(bins)+1));
    [~, largestID] = max(counts);
    idx = find(bins == largestID);

    W_comp = W(idx, idx);
end

function z = eigenvectorCentrality(W)

    [eigenvectors, eigenvalues] = eig(W');
    [~, idx] = max(real(diag(eigenvalues)));

    z = eigenvectors(:, idx);

    if sum(z) < 0
        z = -z;
    end

    z = z / norm(z);
end

function z = katzCentrality(W, beta, mu)
    lambda = eigs(W', 1);


    z = (eye(size(W)) - (1-beta)*W'/lambda)\(beta*mu); %\ is inverse
end

function z = iterativePageRank(W, beta, mu)
    [w_out, ~] = degree(W);
    w_out(w_out == 0) = 1;
    
    D = sparse(diag(w_out));
    P = D \ W;

    mu = mu / sum(mu);
    z = mu; %initial value

    for k = 1:1000
        z_new = (1 - beta) * P' * z + beta * mu;
        
        if norm(z_new - z) < 1e-8
            break;
        end
        
        z = z_new;
    end
    
    z = z_new;
end

function [x, cont] = consensus(W, i1, i2, cont_i)
    [w_out, ~] = degree(W);
    w_out(w_out == 0) = 1;
    
    D = sparse(diag(w_out));
    P = D \ W;
    
    x = ones(length(w_out), 1) * 0.5;
    x(i1) = 1;
    x(i2) = 0;

    cont = zeros(1,10000);

    for k = 1:10000
        x_new = P * x;
        x_new(i1) = 1;
        x_new(i2) = 0;

        cont(k) = x_new(cont_i);
        
        if norm(x_new - x) < 1e-8
            cont = cont(1:k);
            break;
        end
        
        x = x_new;
    end
    
    x = x_new;
end

%%



%% Centrality

% general
% swe2000, idn2000
% find 3 most central sectors

load("IOdownload.mat");
% econ: countries and years
% io: data
% name: sectors

% W
swe2000 = io.swe2000;
idn2000 = io.idn2000;

%% a) in/out

[swe_out, swe_in] = degree(swe2000);
[idn_out, idn_in] = degree(idn2000);

disp("Swe in")
top_swe_in = findLargest(swe_in, 3);
disp(name(top_swe_in(:,1)));
disp(top_swe_in(:,2)');
disp("Swe out")
top_swe_out = findLargest(swe_out, 3);
disp(name(top_swe_out(:,1)));
disp(top_swe_out(:,2)');
disp("Idn in")
top_idn_in = findLargest(idn_in, 3);
disp(name(top_idn_in(:,1)));
disp(top_idn_in(:,2)');
disp("Idn out")
top_idn_out = findLargest(idn_out, 3);
disp(name(top_idn_out(:,1)));
disp(top_idn_out(:,2)');

%% b) eigenvector
swe_comp = largestComponent(swe2000);
idn_comp = largestComponent(idn2000);
swe_eigen = eigenvectorCentrality(swe_comp);
idn_eigen = eigenvectorCentrality(idn_comp);
top_swe_eigen = findLargest(swe_eigen, 3);
top_idn_eigen = findLargest(idn_eigen, 3);
disp("top swe eigen")
disp(name(top_swe_eigen(:,1)));
disp(top_swe_eigen(:,2)');
disp("top idn eigen")
disp(name(top_idn_eigen(:,1)));
disp(top_idn_eigen(:,2)');

%% c) Katz
n = length(swe2000);
beta = 0.15;
mu1 = ones(n, 1);
mu2 = zeros(n, 1);
mu2(31) = 1;

swe_katz_1 = katzCentrality(swe2000, beta, mu1);
swe_katz_2 = katzCentrality(swe2000, beta, mu2);
idn_katz_1 = katzCentrality(idn2000, beta, mu1);
idn_katz_2 = katzCentrality(idn2000, beta, mu2);

top_swe_katz_1 = findLargest(swe_katz_1, 3);
top_swe_katz_2 = findLargest(swe_katz_2, 3);
top_idn_katz_1 = findLargest(idn_katz_1, 3);
top_idn_katz_2 = findLargest(idn_katz_2, 3);

disp("top swe katz ones")
disp(name(top_swe_katz_1(:,1)));
disp(top_swe_katz_1(:,2)');
disp("top swe katz single one")
disp(name(top_swe_katz_2(:,1)));
disp(top_swe_katz_2(:,2)');
disp("top idn katz ones")
disp(name(top_idn_katz_1(:,1)));
disp(top_idn_katz_1(:,2)');
disp("top idn katz single ones")
disp(name(top_idn_katz_2(:,1)));
disp(top_idn_katz_2(:,2)');
%% Influence

users = load('users.mat', '-ascii');
twitter = load('twitter.mat', '-ascii');
n = length(users);
W = spconvert([twitter; n n 0]); %%Would overwrite a selfloop in final node

%% a)
beta = 0.15;
mu = ones(length(W), 1);

pageRank = iterativePageRank(W, 0.15, mu);
top_pageRank = findLargest(pageRank, 5);

disp("top iterative page rank")
disp(top_pageRank(:,1)');
disp(top_pageRank(:,2)');

%% b)

i = [9, 1, 112];

[opinion, cont] = consensus(W, i(1), i(2), i(3));
t = 1:length(cont);
figure
plot(t, cont)
xlabel('Time step')
ylabel(['State of node ', num2str(i(3))])

title(sprintf('Consensus dynamics: node %d over time\nStubborn nodes: %d and %d', ...
    i(3), i(1), i(2)))

%% c)

i = [2, 26, 112];

[opinion, ~] = consensus(W, i(1), i(2), i(3));
t = 1:length(cont);
figure
histogram(opinion)

title(sprintf('Stubborn nodes: %d and %d', ...
    i(1), i(2)))

