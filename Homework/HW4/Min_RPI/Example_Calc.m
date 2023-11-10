% Approximation of Minimal Robust Positively Invariant Set
% Rakovic, Kerrigan, Kouramas, Mayne; TAC, 2005.

A = [1 1; 0 1];
B = [1; 1];

K = -[1.17 1.03];

Ak = A+B*K;

w_max = 1*ones(2,1);
w_min = -w_max;
W_set = Polyhedron('lb',w_min,'ub',w_max);

F = W_set;
F.minHRep;
figure;hold on
F.plot;
drawnow
for i = 1:5
F = plus(F,affineMap(W_set,Ak^i));
F.minHRep;
F.plot; alpha(0.1)
drawnow
end
hold off


%% Approximation of F_inf
% epsilon = 5e-5; % If you use this one, you have to zoom in to see the
% difference in the sets
epsilon = 5e-1;

F_approx = Approx_RPI(Ak,W_set,epsilon);

figure; hold on
F_approx.plot('color', 'blue'); alpha(1)

F.plot; alpha(0.2)
for i = 1:size(F.V,1)
    plot(F.V(i,1),F.V(i,2),'o')
end