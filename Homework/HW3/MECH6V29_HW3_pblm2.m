%% MECH 6V29 - MPC - Homework 3
%% Problem 2

%% 2a
A = [1, 1;
    0, 1];
B = [0.5;
    1];
C = eye(3,2);%<--- [1,0;0,1;0,0]
D(3,1) = 1; %<--- [0;0;1]
sys = ss(A,B,C,D,1)

N = 10;

% sizes
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

% controller
K = -acker(A,B,zeros(nx,1));

% sets
X = Polyhedron('A',[eye(nx);-eye(nx)], 'b', ones(2*nx,1));
U = Polyhedron('A', [1;-1], 'b', ones(2*nu,1));
W = B*Polyhedron('A',[1;-1],'b',[0.3;0.3]);


%% 2b - RPI Set
A_K = A+B*K;
F = W;
F.minHRep;
fig = figure; hold on;
F.plot
% drawnow
for i = 1:5
    F = F + (A_K)^(i)*W;
    F.minHRep;
    F.plot; alpha(0.1);
    % drawnow
end
hold off
saveas(fig,strcat('figs',filesep,'pblm2b_1','.png'));

% Approx code
epsilon = 1; %<== all epsilon smaller then 1 appear to make it good
F_approx = Approx_RPI(A_K,W,epsilon);

% Plot
fig = figure; hold on
F_approx.plot('color', 'blue'); alpha(1)

F.plot; alpha(0.2);
for i = 1:size(F.V,1)
    plot(F.V(i,1),F.V(i,2),'o');
end
saveas(fig,strcat('figs',filesep,'pblm2b_2','.png'))

%% 2c - tightened state/input sets
Z = F_approx;

X_bar = X - Z; X_bar.minHRep;
U_bar = U - K*Z; U_bar.minHRep;

%% 2d - new controller







% Y_{1} = Y - (C+D*K)*W;
% for j = 2:N
%     Y_{j} = Y_{j-1} - (C+D*K)*(A_K)^(j-1)*W;
% end
% 
% for robustFlag = [true, false]
% %% 1d ----- Setup Controller
% P=0;
% Q = 1e-3*eye(nx);
% R = 100;
% 
% yalmip('clear'); clear('controller');
% u_ = sdpvar(repmat(nu,1,N),ones(1,N));
% x_ = sdpvar(repmat(nx,1,N+1),ones(1,N+1));
% 
% constraints = []; objective = 0;
% for k = 1:N
%     objective = objective + x_{k}'*Q*x_{k} + u_{k}'*R*u_{k};
%     constraints = [constraints, x_{k+1} == A*x_{k} + B*u_{k}];
%     if robustFlag
% 
%         constraints = [constraints, Y_{k}.A*(C*x_{k}+D*u_{k})<= Y_{k}.b];
%     else
%         constraints = [constraints, Y.A*(C*x_{k}+D*u_{k})<= Y.b];
%     end
% end
% constraints = [constraints, x_{k+1} == 0];
% objective = objective + x_{k+1}'*P*x_{k+1};
% 
% opts = sdpsettings;
% controller = optimizer(constraints,objective,opts,x_{1},u_{1});
% 
% 
% % simulate and plot
% fig = figure(...
%         WindowStyle="normal",...
%         Position=[0 0 750 500]);
% hold on
% for i = 1:100
%     rng(i);
%     x0 = zeros(nx,1); tf = 100;
%     V = num2cell(0.6*rand(nx,tf)-0.3);
%     [X{i},U{i},~] = run_sim(A,B,V,controller, x0, tf);
%     k_fail = find(~isfinite(U{i}),1,"first");
%     plot(X{i}(1,:),'k')
%     plot(k_fail,X{i}(1,k_fail),'ko')
% end
% yline(1,'k'); yline(-1,'k');
% ylabel('Position');
% xlabel('Time');
% title(sprintf('robustFlag = %d',robustFlag))
% saveas(fig,strcat('figs',filesep,sprintf('pblm1_robust=%d',robustFlag),'.png'));
% 
% end





%% Local functions
% function controller = mpc_yalmip_controller(A,B,P,Q,R,N,cons,cons_f)
%     yalmip('clear')
%     nx = size(A,1);
%     nu = size(B,2);
% 
%     u_ = sdpvar(repmat(nu,1,N),ones(1,N));
%     x_ = sdpvar(repmat(nx,1,N+1),ones(1,N+1));
%     s_ = sdpvar(ones(1,N+1),ones(1,N+1));
% 
%     constraints = [];
%     objective = 0;
%     for k = 1:N
%         objective = objective + x_{k}'*Q*x_{k} + u_{k}'*R*u_{k} + s_{k};
%         constraints = [constraints, s_{k} >= 0];
%         constraints = [constraints, x_{k+1} == A*x_{k} + B*u_{k}];
%         constraints = [constraints, cons(x_{k+1},u_{k},s_{k})];
%     end
%     constraints = [constraints,cons_f(x_{k+1},u_{k},s_{k})];
%     objective = objective + x_{k+1}'*P*x_{k+1};
% 
%     opts = sdpsettings;
%     controller = optimizer(constraints,objective,opts,x_{1},u_{1});
% end

function [X,U,diagnostics_] = run_sim(A,B,V,controller,x0, tf)
    
    X_{tf+1} = []; U_{tf} = []; diagnostics_{tf} = [];
    X_{1} = x0;
    for k = 1:tf
        [U_{k},diagnostics_{k}] = controller{X_{k}};
        X_{k+1} = A*X_{k} + B*U_{k} + B*V{k};
    end
    X = [X_{:}]; U = [U_{:}];
end

% function fig = plot_trajectory(X, U)
%     fig = figure(...
%         WindowStyle="normal",...
%         Position=[0 0 750 500]);
%     hold on; grid on;
%     subplot(2,1,1);
%     stairs(X')
%     title('State Trajectory')
%     legend({'x_1','x_2'})
%     subplot(2,1,2);
%     stairs(U');
%     title('Input Trajectory')
%     legend({'u_1'})
% end

