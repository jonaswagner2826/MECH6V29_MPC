%% System Definition
% System dimensions
nx = 4; % Number of states
nu = 1; % Number of inputs (controllable)
nd = 2; % Number of disturbances
ny = 3; % Number of outputs

% Model parameters
Ks = 900;   % Suspension Stiffness (N/m) 
Kt = 2500;  % Tire stiffness (N/m)
Ms = 2.45;  % Sprung Mass (kg) 
Mu = 1;     % Unsprung Mass (kg)
Bs = 7.5;   % Suspension Inherent Damping coefficient (sec/m)
Bt = 5;    % Tire Inhenrent Damping coefficient (sec/m)

% Continuous-time state-space model
% x1 = z1-z0 (length of tire spring)
% x2 = z1_dot (velocity of unsprung mass)
% x3 = z2-z1 (length of suspension spring)
% x4 = x2_dot (velocity of sprung mass)
% u1 = F (active suspension force)
% d1 = z0_dot (rate of change of road height)
% d2 = z0 (road height)
% y1 = z1 (unsprung mass height)
% y2 = z2 (sprung mass height)
% y3 = z2_dotdot (acceleration of sprung mass)
Ac = [0 1 0 0;...
     -Kt/Mu -(Bs+Bt)/Mu Ks/Mu Bs/Mu;...
      0 -1 0 1;...
      0 Bs/Ms -Ks/Ms -Bs/Ms];
Bc = [0; -1/Mu; 0; 1/Ms];
Vc = [-1 0; Bt/Mu 0; 0 0; 0 0];
C = [1 0 0 0;...
      1 0 1 0;...
      0 Bs/Ms -Ks/Ms -Bs/Ms];
D = [0; 0; 1/Ms];
W = [0 1; 0 1; 0 0];
% LTI system 
sys_c = ss(Ac,[Bc Vc],C,[D W]);

%% Constraints
% Maximum tire deflection +/- 0.01 meters
% Maximum unsprung mass velocity is +/- 1 m/s
% Maximum suspension deflection +/- 0.03 meters
% Maximum unsprung mass velocity is +/- 1 m/s
% Maximum force is +/- 30 N
% Maximum change in road profile velocity +/- 0.5 meters/second
% Maximum change in road height +/- 0.02 meters
bounds.x_ub = [0.01; 1; 0.03; 1];
bounds.x_lb = -bounds.x_ub;
bounds.u_ub = 30;
bounds.u_lb = -bounds.u_ub;
bounds.d_ub = [0.5; 0.02];
bounds.d_lb = -bounds.d_ub;

X_set = Polyhedron('lb',bounds.x_lb,'ub',bounds.x_ub);
U_set = Polyhedron('lb',bounds.u_lb,'ub',bounds.u_ub);
D_set = Polyhedron('lb',bounds.d_lb,'ub',bounds.d_ub);

%% Open-loop Simulation (no active suspension force)
tmax = 5;                               % simulation time length (must be less than 10 - will run out of road)
scale = 1;                              % Scaling of road height
x0 = zeros(nx,1);                       % initial condition
v = 10;                                 % vehicle velocity (m/s) (will change the simulation time step)
% Chose your road profile
% load IRI_737b                         % road profile data (realistic road)
% load roadRamps                        % road profile data (multiple large ramps)
% load roadSingleBump                   % road profile data (one bump)
load roadBumpHole                       % road profile data (one bump and one hole)
dx = road_x(2) - road_x(1);             % spacial step for input data
dt = dx/v;                              % simulation time step
t = 0:dt:tmax;                          % simulation time
z0 = road_z(1:tmax*v/dx+1)*scale;       % road height
z0dot = [0 diff(z0)/dt];            	% road profile velocity
pos = v*t;                              % position along the road
sys_d = c2d(sys_c,dt);                  % Discrete-time state-space model for simulation

% Collect inputs/disturbances
u_OL = zeros(length(t),nu+nd);          % Zero force
u_OL(:,2) = z0dot;                      % Rate of change of road height
u_OL(:,3) = z0;                         % Road height

% Simulate open-loop system
[y_OL,~,x_OL] = lsim(sys_d,u_OL,t,x0);
% Plot all simulation data
plotActiveSuspension(t,u_OL,x_OL,y_OL,bounds)
% Animate simulation data (can comment out if you dont want to animate)
animateActiveSuspension(t,u_OL,y_OL,road_x,road_z,pos,scale)

%% Discrete-time LQR Control (Optional)
% LQR cost function matrices 
Q = ; % DESIGN THIS
R = ; % DESIGN THIS
% LQR design using only the controllable input (force)
[K_LQR,~,~] = dlqr(sys_d.A,sys_d.B(:,1),Q,R);
K_LQR = -K_LQR; % For positive u = Kx
% Closed-loop system under LQR
sys_LQR = ss(sys_d.A+sys_d.B(:,1)*K_LQR,sys_d.B(:,2:end),sys_d.C+sys_d.D(:,1)*K_LQR,sys_d.D(:,2:end),dt);

% Simulate open-loop system
[y_LQR,~,x_LQR] = lsim(sys_LQR,u_OL(:,2:end),t,x0);
u_LQR = u_OL;
u_LQR(:,1) = (K_LQR*x_LQR')';

% Plot all simulation data
plotActiveSuspension(t,u_LQR,x_LQR,y_LQR,bounds)
% Animate simulation data (can comment out if you dont want to animate)
animateActiveSuspension(t,u_LQR,y_LQR,road_x,road_z,pos,scale)

%% Reachability Analysis (Optional)

%% Invariant Set (Optional)

%% MPC Controller Design
% Discrete-time step for MPC (must be an integer multiple of dt)
dt_MPC = ; % DESIGN THIS
% Resample discrete-time model with MPC time step
sys_MPC = d2d(sys_d, dt_MPC);

% Prediction horizon
N = ; % DESIGN THIS
% Optimization variables (Choose your Yalmip sdpvar varaibles)

constraints = [];
objective = 0;
for k = 1:N
    % DESIGN YOUR CONTROLLER OBJECTIVE FUNCTION AND CONSTRAINTS
end

controller = optimizer(constraints,objective,sdpsettings('solver','gurobi'),{},[]);
%% Closed-loop MPC Simulation
x_sim = [x0];           % Store states
u_sim = [];             % Store inputs
y_sim = [];             % Store outputs
s_sim = [];             % Store slack (if you use it)
errors = [];            % Store diagnostics
tCalcs = [];            % Store MPC computation times
z0dot2 = [z0dot z0dot]; % Extend for prediction horizon
% Choose what the MPC controller knows (0 - no disturbance, 1 - current
% distrubance, 2 - full disturbance preview)
dFlag = 0;
for i = 1:length(t)-1
    % Call controller every dt_MPC seconds
    if mod(i-1,dt_MPC/dt) < eps 
        if dFlag == 0
            d_MPC = mat2cell(zeros(1,N),1,ones(1,N));
        elseif dFlag == 1
            d_MPC = mat2cell(repmat(z0dot2(i+1),1,N),1,ones(1,N));
        elseif dFlag == 2
            d_MPC = mat2cell(z0dot2(i+1:dt_MPC/dt:i+dt_MPC/dt*N),1,ones(1,N));
        end
        % Modify this section as needed depending on your inputs/outputs
        % Define controller inputs
        
        tic
        [out,diagnostics] = controller{};
        tCalc = toc;
        % Unpack controller outputs with optimal input at first time step
        % denoted by the variable u (to be used below)
        
        tCalcs = [tCalcs tCalc];
        errors = [errors, diagnostics];
    end
    % Update model every dt seconds
    u_sim = [u_sim u];
    y_sim = [y_sim sys_d.C*x_sim(:,end) + sys_d.D(:,1)*u + sys_d.D(:,2:3)*[z0dot(i+1); z0(i+1)]];
    x_sim = [x_sim sys_d.A*x_sim(:,end) + sys_d.B(:,1)*u + sys_d.B(:,2:3)*[z0dot(i+1); z0(i+1)]];
end
% Display diagnostics (0 = no errors)
max(errors)
% Rotate for plotting
x_sim = x_sim';
u_sim = u_sim';
u_sim = [u_sim;u_sim(end)];
u_sim = [u_sim z0dot' z0'];
y_sim = y_sim';
y_sim = [y_sim;y_sim(end,:)];

% Plot all simulation data
plotActiveSuspension(t,u_sim,x_sim,y_sim,bounds)
% Animate simulation data (can comment out if you dont want to animate)
animateActiveSuspension(t,u_sim,y_sim,road_x,road_z,pos,scale)
