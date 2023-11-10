% Approximation of Minimal Robust Positively Invariant Set
% Rakovic, Kerrigan, Kouramas, Mayne; TAC, 2005.

A = [1 1; 0 1];
B = [1; 1];

K = -[1.17 1.03];

Ak = A+B*K;

w_max = 1*ones(2,1);
w_min = -w_max;
W_set.c = zeros(2,1);
W_set.G = diag(w_max);

F = W_set;
figure;hold on

Box = Polyhedron('lb',-ones(size(F.G,2),1),'ub',ones(size(F.G,2),1));
F_Hrep = plus(F.c,affineMap(Box,F.G));
F_Hrep.plot;
drawnow
for i = 1:5
    F.c = F.c + Ak^i*W_set.c;
    F.G = [F.G Ak^i*W_set.G];
    Box = Polyhedron('lb',-ones(size(F.G,2),1),'ub',ones(size(F.G,2),1));
    F_Hrep = plus(F.c,affineMap(Box,F.G));
    F_Hrep.plot; alpha(0.1)
    drawnow
end
hold off


%% Approximation of F_inf
% epsilon = 5e-5; % If you use this one, you have to zoom in to see the
% difference in the sets
epsilon = 5e-1;

F_approx = Approx_RPI(Ak,W_set,epsilon);

figure;hold on
Box = Polyhedron('lb',-ones(size(F_approx.G,2),1),'ub',ones(size(F_approx.G,2),1));
F_approx_Hrep = plus(F_approx.c,affineMap(Box,F_approx.G));
F_approx_Hrep.plot('color', 'blue'); alpha(1)

F_Hrep.plot; alpha(0.2)
for i = 1:size(F_Hrep.V,1)
    plot(F_Hrep.V(i,1),F_Hrep.V(i,2),'o')
end