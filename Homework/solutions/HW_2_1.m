%% Define System Matrices and discrete-time system
n = 2;
m = 1;
A = [4/3 -2/3; 1 0];
B = [1; 0];
C = [-2/3 1];
D = [0];
x0 = [3; 3];           % Initial condition

%% Define Objective Function Parameters
Q = eye(n);             % State penalties
R = eye(m);             % Input penalties
P = 0;                  % Terminal cost
N = 10;                 % Prediction horizon

%% Controller Formulation
u = sdpvar(repmat(m,1,N),repmat(1,1,N));
x = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
    constraints = [constraints, C*x{k+1} <= 0.5];
end
constraints = [constraints, x{N+1} == 0];

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x{1},u{1});

%%
k_sim = 10;
x_sim = [x0];
u_sim = [];
y_sim = [];
errors = [];
for i = 1:k_sim
    [u,diagnostics] = controller{x_sim(:,i)};
    errors(i) = diagnostics;
    u_sim = [u_sim u];
    y_sim = [y_sim C*x_sim(:,i)];
    x_sim = [x_sim A*x_sim(:,end) + B*u_sim(:,end)];
end
errors

figure;
subplot(3,1,1); hold on
stairs(0:k_sim,x_sim(1,:))
stairs(0:k_sim,x_sim(2,:))
xlabel('Time (s)')
ylabel('States')
legend('x_1','x_2')
% ylim([-6 10])
subplot(3,1,2);
stairs(0:k_sim-1,y_sim)
xlabel('Time (s)')
ylabel('Output')
% ylim([-10 5])
subplot(3,1,3);
stairs(0:k_sim-1,u_sim)
xlabel('Time (s)')
ylabel('Input')
% ylim([-10 5])


%% Controller Formulation
u = sdpvar(repmat(m,1,N),repmat(1,1,N));
x = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));
s = sdpvar(1,1);

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
    constraints = [constraints, C*x{k+1} <= 0.5 + s;];
end
constraints = [constraints, x{N+1} == 0];
constraints = [constraints, s >= 0];
objective = objective + 1e6*s^2;

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x{1},u{1});

%%
k_sim = 10;
x_sim = [x0];
u_sim = [];
y_sim = [];
errors = [];
for i = 1:k_sim
    [u,diagnostics] = controller{x_sim(:,i)};
    errors(i) = diagnostics;
    u_sim = [u_sim u];
    y_sim = [y_sim C*x_sim(:,i)];
    x_sim = [x_sim A*x_sim(:,end) + B*u_sim(:,end)];
end
errors

figure;
subplot(3,1,1); hold on
stairs(0:k_sim,x_sim(1,:))
stairs(0:k_sim,x_sim(2,:))
xlabel('Time (s)')
ylabel('States')
legend('x_1','x_2')
% ylim([-6 10])
subplot(3,1,2);
stairs(0:k_sim-1,y_sim)
xlabel('Time (s)')
ylabel('Output')
% ylim([-10 5])
subplot(3,1,3);
stairs(0:k_sim-1,u_sim)
xlabel('Time (s)')
ylabel('Input')
% ylim([-10 5])



