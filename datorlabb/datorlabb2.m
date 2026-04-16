%% Pre-processor

geom2;

%A = 10;
k = 1;
%L = 6;
%Q = 100;
%nen = 1;
%nelm = 100; %edit this

%Le = 6/nelm;
%ndof = nelm+1;

%Coord = (2:Le:8)';
%Coord = [2; 4; 6; 8];

%Edof = zeros(nelm, 3);
%Edof(:,1) = 1:nelm;
%Edof(:,2) = 1:ndof-1;
%Edof(:,3) = 2:ndof;
%Edof = [1 1 2; 2 2 3; 3 3 4];

%Dof = (1:ndof)';
%Dof = [1; 2; 3; 4];

%[Ex] = coordxtr(Edof, Coord, Dof, nen);
%Ex = Coord;

K = zeros(ndof);
F = zeros(ndof, 1);

%F(ndof) = F(ndof) + 15*A;
%bc = [1 0];

%% Solver

for elnr = 1:nelm
    Ke = flw2te(ex(elnr,:), ey(elnr,:), 1, [k 0; 0 k]);
    K = assem(edof(elnr,:), K, Ke);
end

a = solve(K, F, bc);

%% Post-processor


%xplot = 1:ndof;
%plot(xplot, a);
ed = extract(edof, a);
es = zeros(nelm, 2);
et = zeros(nelm, 2);
figure
hold on
for elnr = 1:nelm
    [es(elnr,:), et(elnr,:)] = flw2ts(ex(elnr,:), ey(elnr,:), [k 0; 0 k], ed(elnr,:));
end
es = es/norm(es)*2;
es(es < 0) = 0;
for elnr = 1:nelm
    fill(ex(elnr,:), ey(elnr,:), [es(elnr, 1) es(elnr, 2) 0]);
end
hold off

%fill(ex(1:2,:), ey(1:2,:), [1 0 0; 0 1 0])
%fill(ex, ey, [1 0 0])
