%% MPC HW 1 - Problem 1
% clear
% close all
subfolder = fileparts(mfilename('fullpath'));
if ~isfolder(strcat(subfolder,filesep,'figs'))
    mkdir(strcat(subfolder,filesep,'figs'));
end

% Problem Information
A = [4/3, -2/3; 1, 0];
B = [1; 0];
C = [-2/3, 1];
D = 0;
dt = 1;
sys = ss(A,B,C,D,dt);

% Size parameters
nx = size(A,1);
nu = size(B,2);

% MPC Parameters
N = 5;
Q = eye(nx);
R = eye(nu);
P = 0;

%% Part 1a
% Controller setup
cons = @(x,u) [];
controller = mpc_yalmip_controller(A, B, P, Q, R, N, cons);

% Simulation
x0 = [3;3]; tf = 10;
[X, U] = run_sim(A,B,controller, x0, tf);
Y = C*X;

% Results
figName = 'pblm1a';
fig = plot_trajectory(X, Y, U);
saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

maxOutput = max(Y); 
fprintf('Max Output: %f\n',maxOutput);


%% Part 1b
% Controller setup
cons = @(x,u) (C*x)<=0.5;
R = 0;
controller = mpc_yalmip_controller(A, B, P, Q, R, N, cons);

% Simulation
x0 = [3;3]; tf = 10;
[X, U] = run_sim(A,B,controller, x0, tf);
Y = C*X;

% Results
figName = 'pblm1b';
fig = plot_trajectory(X, Y, U);
subplot(3,1,2); yline(0.5,'--','Output Constraint')
saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')



















% [y,t,x] = step(sys);

% % Output Response
% figName = 'pblm1a_fig1';
% fig = figure(WindowStyle="normal");
% hold on; grid on;
% stairs(t,y,"DisplayName",'Output');
% title('Output Response')
% xlabel('Time (s)')
% ylabel('Output')
% saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

% % State Response
% figName = 'pblm1a_fig2';
% fig = figure(WindowStyle="normal");
% hold on; grid on;
% stairs(t,x(:,1),"DisplayName",'x_1');
% stairs(t,x(:,2),"DisplayName",'x_2');
% legend
% title('State Response')
% xlabel('Time (s)')
% ylabel('')
% saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

% % State Response
% figName = 'pblm1a_fig3';
% fig = figure(WindowStyle="normal");
% hold on; grid on;
% stairs(x(:,1),x(:,2))
% title('State Traces')
% xlabel('x_1')
% ylabel('x_2')
% saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

% %% Part 1b
% % MPC Parameters (standard optimization)
% Q = C'*C + 1e-3*eye(nx);
% R = 1e-3;
% P = Q;
% N = 5; % prediction horizon


% %% Batch method
% % Construct S_x,and S_u
% for i = 1:N+1
%     S_x_{i} = A^(i-1);
%     for j = 1:nu
%         S_u__{j} = eye(nu)*A^(i-2+j-1)*B;
%     end
%     S_u_{i} = horzcat(S_u__{:});
% end
% S_x = vertcat(S_x_{:});
% S_u = tril(vertcat(S_u_{:}),-1);

% % Construct Qbar and Rbar
% Qbar = blkdiag(kron(Q,eye(N)),P);
% Rbar = kron(R,eye(N));

% % H, F, Y
% H = S_u'*Qbar*S_u + R; 
% F = S_x'*Qbar*S_u;
% Y = S_x'*Qbar*S_x;

% % Results
% Ustar = @(x_0) -H\F'*x_0;
% ustar = @(x_0) [eye(nx),zeros(nx,nx*(N-1))]*-H\F'*x_0;
% K_0_batch = -H\F';

% disp('Batch Method:')
% disp('K = '); disp(K_0_batch);

% %% Dynamic Programing Method
% F_fun = @(P_k) -(B'*P_k*B+R)\B'*P_k*A;
% P_fun = @(P_kp1) A'*P_kp1*A + Q - A'*P_kp1*B*inv(B'*P_kp1*B+R)*B'*P_kp1*A;

% P_{N} = Q;

% for k = N-1:-1:1
%     P_{k} = P_fun(P_{k+1});
%     F_{k} = F_fun(P_{k});
% end
% P_0 = P_fun(P_{1});
% F_0 = F_fun(P_0);

% ustar = @(x) F_0*x;
% K_0_dynprog = F_0;

% disp('Dynamic Programing Method:')
% disp('K = '); disp(K_0_dynprog);

% % results
% disp('They are not the same, or at least not with a time-horrizon on $N=5$')

% %% Part c
% disp('Batch version')
% A_K_batch = A+B*K_0_batch
% eig_batch = eig(A_K_batch)

% disp('The system is not closed-loop stable acording to this as there is an eigen value >= 1')

% disp('Dynamic Programing Method')
% A_K_dynprog = A+B*K_0_dynprog
% eigh_dynprog = eig(A_K_dynprog)

% disp('The system is not closed-loop stable acording to this as there is an eigen value >= 1')

% %% Part d
% K_lqr = -dlqr(A,B,Q,R)
% A_K_lqr = A+B*K_lqr
% eig_lqr = eig(A_K_lqr)

% disp('LQR is stable')

% %% Part e

% for N = 1:20
%     F_fun = @(P_k) -(B'*P_k*B+R)\B'*P_k*A;
%     P_fun = @(P_kp1) A'*P_kp1*A + Q - A'*P_kp1*B*...
%         inv(B'*P_kp1*B+R)*B'*P_kp1*A;
    
%     P_{N} = Q;    
%     for k = N-1:-1:1
%         P_{k} = P_fun(P_{k+1});
%         F_{k} = F_fun(P_{k});
%     end
%     P_0 = P_fun(P_{1});
%     F_0 = F_fun(P_0);

%     K_0_dynprog_{N} = F_0;
%     A_K_dynprog_{N} = A+B*K_0_dynprog_{N};
%     eig_dynprog_{N} = eig(A+B*K_0_dynprog_{N});
%     max_eig_{N} = max(abs(eig_dynprog_{N}));
% end

% max_eig = [max_eig_{:}];

% figName = 'pblm1e';
% fig = figure(WindowStyle="normal");
% hold on; grid on;
% plot(1:N,max_eig)
% yline(max(eig_lqr))
% legend({'Dyn Prog', 'LQR'})
% title('Stability Comparrision')
% xlabel('Prediction Horrizon (N)')
% ylabel('max(abs(\lamba))')
% saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')



% End
close all



%% Local Functions
function controller = mpc_yalmip_controller(A,B,P,Q,R,N, cons)
    yalmip('clear')
    nx = size(A,1);
    nu = size(B,2);
    
    u_ = sdpvar(repmat(nu,1,N),ones(1,N));
    x_ = sdpvar(repmat(nx,1,N+1),ones(1,N+1));
    
    constraints = [];
    objective = 0;
    for k = 1:N
        objective = objective + x_{k}'*Q*x_{k} + u_{k}'*R*u_{k}; %norm(Q*x_{k},2) + norm(R*u_{k},2);
        constraints = [constraints, x_{k+1} == A*x_{k} + B*u_{k}];
        constraints = [constraints, cons(x_{k+1},u_{k})];
        % constraints = [constraints, x_con(x_{k})];
    end
    objective = objective + x_{k+1}'*P*x_{k+1};
    
    opts = sdpsettings;
    controller = optimizer(constraints, objective,opts,x_{1},u_{1});
end

function [X,U] = run_sim(A,B,controller,x0, tf)
    
    X_{tf+1} = []; U_{tf} = [];
    X_{1} = x0;
    for k = 1:tf
        U_{k} = controller{X_{k}};
        X_{k+1} = A*X_{k} + B*U_{k};
    end
    X = [X_{:}]; U = [U_{:}];
end

function fig = plot_trajectory(X, Y, U)
    fig = figure(WindowStyle="normal");
    hold on; grid on;
    subplot(3,1,1);
    stairs(X')
    title('State Trajectory')
    legend({'x_1','x_2'})
    subplot(3,1,2);
    stairs(Y');
    title('Output Trajectory')
    legend({'y_1'})
    subplot(3,1,3);
    stairs(U');
    title('Input Trajectory')
    legend({'u_1'})
end