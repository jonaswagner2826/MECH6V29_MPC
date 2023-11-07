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
% Output set
Y0 = Polyhedron('lb',-1*ones(n+m,1),'ub',ones(n+m,1));

% Disturbance Set
z_max = 0.3;
Z = Polyhedron('lb',-z_max,'ub',z_max);
W = B*Z;  

Y{1} = Y0;
L{1} = eye(n);
for j = 1:N-1
    Y{j+1} = Y{j} - (C+D*K)*L{j}*W;
    L{j+1} = (A+B*K)*L{j};
end
%% Controller formulation
u_ = sdpvar(repmat(m,1,N),repmat(1,1,N));
x_ = sdpvar(repmat(n,1,N+1),repmat(1,1,N+1));

robustFlag = 1;

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + 1e-3*x_{k}'*x_{k} + 1e2*u_{k}'*u_{k};
    constraints = [constraints, x_{k+1} == A*x_{k} + B*u_{k}];
    if robustFlag == 0
        constraints = [constraints, Y0.H(:,1:end-1)*(C*x_{k}+D*u_{k}) <= Y0.H(:,end)];
    else
        constraints = [constraints, Y{k}.H(:,1:end-1)*(C*x_{k}+D*u_{k}) <= Y{k}.H(:,end)];
    end
end
constraints = [constraints, x_{N+1} == 0];

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),x_{1},u_{1});

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
        [u,diagnostics] = controller{x_sim(:,i)};
        if diagnostics == 0
            u_sim = [u_sim u];
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

mean(Costs) % 257.5789
max(Costs)  % 460.0397