%% Pre-processor

A = 10;
k = 5;
L = 6;
Q = 100;
nen = 1;
nelm = 100; %edit this

Le = 6/nelm;
ndof = nelm+1;

Coord = (2:Le:8)';
%Coord = [2; 4; 6; 8];

Edof = zeros(nelm, 3);
Edof(:,1) = 1:nelm;
Edof(:,2) = 1:ndof-1;
Edof(:,3) = 2:ndof;
%Edof = [1 1 2; 2 2 3; 3 3 4];

Dof = (1:ndof)';
%Dof = [1; 2; 3; 4];

%[Ex] = coordxtr(Edof, Coord, Dof, nen);
Ex = Coord;

K = zeros(ndof);
F = zeros(ndof, 1);

F(ndof) = F(ndof) + 15*A;
bc = [1 0];

%% Solver

for elnr = 1:nelm
    Ke = spring1e(k*A/Le);
    Fe = Q*A*Le/2 * [1; 1];
    [K, F] = assem(Edof(elnr,:), K, Ke, F, Fe);
end

a = solve(K, F, bc);

%% Post-processor


xplot = 1:ndof;
plot(xplot, a);
%ed = extract(Edof, a);
%es = spring1s(k*A/Le, ed);
%eldraw
%eldisp