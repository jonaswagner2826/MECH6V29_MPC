
% syms s t

K = 1e6;
m = 1e-4;
k = 100;
c = 0.1;

% G = (K*m)/(m*s^2 + c*s + k)
% 
% a_z1 = 5*sin(100*t + 0.3);
% a_z2 = 20*sin(1500*t + 0.2);
% 
% A_z1 = fourier(a_z1)
% 
% Y_1 = A_z1*G
% y_z1 = simplify(ilaplace(Y_1))


G = tf(K*m,[m c k])
% t = 0:1e-5:0.5;
% a_z1 = 5*sin(100*t + 0.3);

% y = lsim(G,a_z1,t);

% figure
% hold on
% plot(t,a_z1)
% plot(t,y)

bode(G)