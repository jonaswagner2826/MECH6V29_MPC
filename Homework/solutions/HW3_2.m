%% Define System Matrices and discrete-time system
n = 2;
m = 1;
A = [1 1; 0 1];
B = [0.5; 1];
C = [1 0; 0 1; 0 0];
D = [0; 0; 1];
x0 = zeros(n,1);           % Initial condition
N = 10;                    % Prediction horizon

%% Candidate feedback controller
K = -acker(A,B,zeros(n,1)); % Remember to use negative sign

%% Sets
% State and Input set
X0 = Polyhedron('lb',-1*ones(n,1),'ub',ones(n,1));
U0 = Polyhedron('lb',-1*ones(m,1),'ub',ones(m,1));

% Disturbance Set
z_max = 0.3;
Z = Polyhedron('lb',-z_max,'ub',z_max);
W = B*Z;  

% RPI set
epsilon = 1e0;
E = Approx_RPI(A+B*K,W,epsilon);
% E.volume
figure;plot(E)

X = X0-E;
U = U0-K*E;

%% Controller formulation
u_ = sdpvar(repmat(m,1,N),repmat(1,1,N));
x_ = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));
x0_ = sdpvar(repmat(n,1,1),repmat(1,1,1));

ICFlag = 1;

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + 1e-3*x_{k}'*x_{k} + 1e2*u_{k}'*u_{k};
    constraints = [constraints, x_{k+1} == A*x_{k} + B*u_{k}];
    constraints = [constraints, X.H(:,1:end-1)*x_{k} <= X.H(:,end)];
    constraints = [constraints, U.H(:,1:end-1)*u_{k} <= U.H(:,end)];
end
constraints = [constraints, x_{N+1} == 0];
if ICFlag == 0
    constraints = [constraints, x0_ == x_{1}];
else
    constraints = [constraints, E.H(:,1:end-1)*(x0_ - x_{1}) <= E.H(:,end)];
end

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x0_,[u_{1};x_{1}]);

%% Simulation Results
rng(1)      % Results might vary based on seed
nSim = 100; % Number of simulations
tSim = 100; % Time of simulations in steps

figure;hold on 
figX = gcf;
title('State 1')
plot([0 tSim],[1 1],'--k')
plot([0 tSim],[-1 -1],'--k')

figure;hold on 
figU = gcf;
title('Input')
plot([0 tSim],[1 1],'--k')
plot([0 tSim],[-1 -1],'--k')

Costs = inf(nSim,1);
for j = 1:nSim
    x_sim = [x0];
    u_sim = [];
    z_sim = [];
    Cost = 0;
    for i = 1:tSim
        [sol,diagnostics] = controller{x_sim(:,i)};
        u_star = sol(1:m);
        x_star = sol(m+1:end);
        if diagnostics == 0
            u_sim = [u_sim u_star+K*(x_sim(:,end)-x_star)];
            z = -z_max + 2*z_max*rand;
            z_sim = [z_sim z];
            x_sim = [x_sim A*x_sim(:,end) + B*u_sim(:,end) + B*z];
            Cost = Cost + 1e-3*x_sim(:,end-1)'*x_sim(:,end-1) + 1e2*u_sim(:,end)'*u_sim(:,end);
        else
            i = i - 1;
            Cost = inf;
            break
        end
    end
    Costs(j) = Cost;
    figure(figX)
    plot(0:i,x_sim(1,1:i+1));
    plot(i,x_sim(1,i+1),'o')
    figure(figU)
    plot(0:i-1,u_sim(1,1:i));
    plot(i-1,u_sim(1,i),'o')
    drawnow
end

mean(Costs) % 1.4625e+03
max(Costs)  % 1.9431e+03


