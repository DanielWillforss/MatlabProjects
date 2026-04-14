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
    G = graph(W);
    bins = conncomp(G);

    counts = histcounts(bins, 1:(max(bins)+1));
    [~, largestID] = max(counts);
    idx = find(bins == largestID);

    W_comp = W(idx, idx);
end

function z = eigenvectorCentrality(W)

    [eigenvectors, eigenvalues] = eig(W);
    [~, idx] = max(abs(diag(eigenvalues)));

    z = eigenvectors(:, idx);
    z = z / norm(z);
end

%% Centrality

A = randi(100, 500, 500);   % 500x500 matrix with integers 1–100

[w_out, w_in] = degree(A);

out = findLargest(w_out, 10);
in = findLargest(w_in, 10);

z = eigenvectorCentrality(largestComponent(A));
center = findLargest(z, 10);

%% Influence




