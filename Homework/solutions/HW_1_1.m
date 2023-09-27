%% Define System Matrices and discrete-time system
n = 2;
m = 1;
A = [ 4/3 -2/3; 1 0];
B = [1; 0];
C = [-2/3 1];
D = [0];
dt = 1; 
sys = ss(A,B,C,D,dt);

%% Open-loop step input
[y,t,x] = step(sys);
% Output as a function of time
figure; stairs(t,y);
xlabel('Time (seconds)')
ylabel('Output')
% States as a function of time
figure; hold on
stairs(t,x(:,1));
stairs(t,x(:,2));
xlabel('Time (seconds)')
ylabel('States')
legend('x_1','x_2')
% States
figure;
stairs(x(:,1),x(:,2))
xlabel('x_1')
ylabel('x_2')

%% Define Cost Function Matrices
Q = C'*C + 1e-3*eye(n); % State penalties
R = 1e-3;               % Input penalties
P = Q;                  % Terminal cost
N = 5;                  % Prediction horizon

%% Batch Approach
Sx = eye(n);
Su = zeros(n,m*N);
Q_bar = [];
R_bar = [];
for i = 1:N
    Sx = [Sx; A*Sx(end-n+1:end,:)];
    Su = [Su; A^(i-1)*B Su(end-n+1:end,1:end-m)];
    Q_bar = blkdiag(Q_bar,Q);
    R_bar = blkdiag(R_bar,R);
end
Q_bar = blkdiag(Q_bar,P);

H = Su'*Q_bar*Su + R_bar;
F = Sx'*Q_bar*Su;
Y = Sx'*Q_bar*Sx;

K_batch = -[1 zeros(1,N-1)]*(H\F')


%% Recursive Approach
P_k = P;
for k = 1:N-1
    P_k = A'*P_k*A + Q - A'*P_k*B*inv(B'*P_k*B+R)*B'*P_k*A;
end
K_rec = -inv(B'*P_k*B+R)*B'*P_k*A

%% Closed-loop eigenvalues
K0 = K_rec;
AK0 = A + B*K0;
eig(AK0)

%% Infinite horizon
[K,~,~] = dlqr(A,B,Q,R);
eig(A-B*K)

%% Effect of N
eigData = [];
for N = 1:20
    P_k = P;
    for k = 1:N-1
        P_k = A'*P_k*A + Q - A'*P_k*B*inv(B'*P_k*B+R)*B'*P_k*A;
    end
    K_rec = -inv(B'*P_k*B+R)*B'*P_k*A;
    eigData(N) = max(abs(eig(A + B*K_rec)));
end
figure;hold on
plot(1:N,eigData,'o')
plot([1 N], [max(abs(eig(A-B*K))) max(abs(eig(A-B*K)))],'--k')
xlabel('N')
ylabel('Max. eig.')
%%
% Ac = [0 1; 0 0];
% Bc = [0; 1];
% Cc = [1 0];
% Dc = [0];
% sys_c = ss(Ac,Bc,Cc,Dc);
% dt = 1; 
% sys = c2d(sys_c,dt);
