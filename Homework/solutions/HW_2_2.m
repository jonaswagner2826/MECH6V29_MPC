%% Define System Matrices and discrete-time system
n = 2;
m = 2;
A = [2 1; 0 2];
B = [1 0; 0 1];

x0 = [-1; 0.5];           % Initial condition

%% Define Objective Function Parameters
alpha = 1e-1;
Q = alpha*eye(n);             % State penalties
R = eye(m);             % Input penalties
P = 0;                  % Terminal cost
N = 3;                 % Prediction horizon

%% Controller Formulation
u = sdpvar(repmat(m,1,N),repmat(1,1,N));
x = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
    constraints = [constraints, -1 <= u{k} <= 1];
    constraints = [constraints, [1 0]*x{k+1} <= 5];
end
constraints = [constraints, x{N+1} == 0];

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x{1},u{1});

%%
k_sim = 10;
x_sim = [x0];
u_sim = [];
errors = [];
for i = 1:k_sim
    [u,diagnostics] = controller{x_sim(:,i)};
    errors(i) = diagnostics;
    u_sim = [u_sim u];
    x_sim = [x_sim A*x_sim(:,end) + B*u_sim(:,end)];
end
errors

figure;
subplot(2,1,1); hold on
stairs(0:k_sim,x_sim(1,:))
stairs(0:k_sim,x_sim(2,:))
xlabel('Time (s)')
ylabel('States')
legend('x_1','x_2')
% ylim([-6 10])
subplot(2,1,2); hold on
stairs(0:k_sim-1,u_sim(1,:))
stairs(0:k_sim-1,u_sim(2,:))
xlabel('Time (s)')
ylabel('Input')
legend('u_1','u_2')
% ylim([-10 5])

%% Region of Attraction (Batch)
Ax = [1 0];
bx = [5];
Af = [eye(n);-eye(n)];
bf = zeros(2*n,1);
Au = [eye(m);-eye(m)];
bu = ones(2*m,1);

tic;
Sx = eye(n);
Su = zeros(n,m*N);
Ax_bar = [];
Au_bar = [];
for i = 1:N
    Sx = [Sx; A*Sx(end-n+1:end,:)];
    Su = [Su; A^(i-1)*B Su(end-n+1:end,1:end-m)];
    Ax_bar = blkdiag(Ax_bar,Ax);
    Au_bar = blkdiag(Au_bar,Au);
end
Ax_bar = blkdiag(Ax_bar,Af);
A_bar = blkdiag(Ax_bar,Au_bar);
b_bar = [repmat(bx,N,1);bf;repmat(bu,N,1)];

H = A_bar*[Sx Su; zeros(m*N,n) eye(m*N)];
f = b_bar;
P = Polyhedron('H',[H f]);

X0 = projection(P,1:n);
toc

figure;plot(X0)
xlim([-2 2])
ylim([-2 2])

%% Region of Attraction (Recursive)
tic
AK = Af;
bK = bf;
for i = 1:N
    H = [AK*A AK*B; zeros(2*m,n) Au];
    f = [bK;bu];
    P = Polyhedron('H',[H f]);
    Pre = projection(P,1:n);
    AK = [Pre.H(:,1:end-1); Ax];
    bK = [Pre.H(:,end); bx];
end

X0 = Polyhedron('H',[AK bK]);
toc
figure;plot(X0)
xlim([-2 2])
ylim([-2 2])

