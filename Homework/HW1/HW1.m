%% MPC HW 1
clear
close all

subfolder = fileparts(mfilename('fullpath'));
if ~isfolder('figs'); mkdir('figs'); end



%% Problem 1
A = [4/3, -2/3; 1, 0];
B = [1; 0];
C = [-2/3, 1];
D = 0;
dt = 1;
sys = ss(A,B,C,D,dt);

% MPC Parameters
Q = C'*C + 1e-3*eye(size(C'*C));
R = 1e-3;

%% Part 1a
% Step Response
[y,t,x] = step(sys);

% Output Response
figName = 'pblm1a_fig1';
fig = figure(WindowStyle="normal");
hold on; grid on;
stairs(t,y,"DisplayName",'Output');
title('Output Response')
xlabel('Time (s)')
ylabel('Output')
saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

% State Response
figName = 'pblm1a_fig2';
fig = figure(WindowStyle="normal");
hold on; grid on;
stairs(t,x(:,1),"DisplayName",'x_1');
stairs(t,x(:,2),"DisplayName",'x_2');
legend
title('State Response')
xlabel('Time (s)')
ylabel('')
saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

% State Response
figName = 'pblm1a_fig3';
fig = figure(WindowStyle="normal");
hold on; grid on;
stairs(x(:,1),x(:,2))
title('State Traces')
xlabel('x_1')
ylabel('x_2')
saveas(fig,[subfolder,filesep,'figs',filesep,figName],'png')

%% Part 1b
