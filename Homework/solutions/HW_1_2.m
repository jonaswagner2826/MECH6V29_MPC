%% Define System Matrices and discrete-time system
n = 2;
m = 1;
% A = [1 1; 0 1];
% B = [0; 1];
% C = [1 0];
% D = [0];
A = [4/3, -2/3; 1, 0];
B = [1; 0];
C = [-2/3, 1];
D = 0;
dt = 1; 
sys = ss(A,B,C,D,dt);
x0 = [3; 3];           % Initial condition

%% Define Objective Function Parameters
Q = eye(n);             % State penalties
R = eye(1);%1e-1;               % Input penalties
P = 0;                  % Terminal cost
N = 5;                  % Prediction horizon

%% Controller Formulation
u = sdpvar(repmat(m,1,N),repmat(1,1,N));
x = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
    % constraints = [constraints, -1 <= u{k} <= 1];
    % constraints = [constraints, -1 <= x{k}(2) <= 1];
end
objective = objective + x{k+1}'*P*x{k+1};

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x{1},u{1});

%%
x_sim = [x0];
u_sim = [];
for i = 1:10
    u_sim = [u_sim controller{x_sim(:,i)}];
    x_sim = [x_sim A*x_sim(:,end) + B*u_sim(:,end)];
end

figure;
subplot(2,1,1); hold on
stairs(0:10,x_sim(1,:))
stairs(0:10,x_sim(2,:))
xlabel('Time (s)')
ylabel('States')
legend('x_1','x_2')
ylim([-6 10])
subplot(2,1,2);
stairs(0:9,u_sim)
xlabel('Time (s)')
ylabel('Input')
ylim([-10 5])


