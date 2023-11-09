function plotActiveSuspension(t,u,x,y,bounds)

set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 14)
set(0,'defaultTextFontSize' , 14)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

% figure('WindowStyle','normal','Position',[0,0,1000,1200])
% figure(1,'WindowStyle','normal','Position',[0,0,1000,1200])
figure(1); clf;
subplot(4,2,1);hold on
plot(t,x(:,1),'k-')
plot([t(1) t(end)],[bounds.x_ub(1) bounds.x_ub(1)],'k--')
plot([t(1) t(end)],[bounds.x_lb(1) bounds.x_lb(1)],'k--')
xlim([0 5])
title('$x_1$ - tire deflection')
subplot(4,2,2);hold on
plot(t,x(:,2),'k-')
plot([t(1) t(end)],[bounds.x_ub(2) bounds.x_ub(2)],'k--')
plot([t(1) t(end)],[bounds.x_lb(2) bounds.x_lb(2)],'k--')
title('$x_2$ - unsprung mass velocity')
xlim([0 5])
subplot(4,2,3);hold on
plot(t,x(:,3),'k-')
plot([t(1) t(end)],[bounds.x_ub(3) bounds.x_ub(3)],'k--')
plot([t(1) t(end)],[bounds.x_lb(3) bounds.x_lb(3)],'k--')
xlim([0 5])
title('$x_3$ - suspension deflection')
subplot(4,2,4);hold on
plot(t,x(:,4),'k-')
plot([t(1) t(end)],[bounds.x_ub(4) bounds.x_ub(4)],'k--')
plot([t(1) t(end)],[bounds.x_lb(4) bounds.x_lb(4)],'k--')
title('$x_4$ - sprung mass velocity')
xlim([0 5])
subplot(4,2,5);hold on
plot(t,u(:,1),'k-')
plot([t(1) t(end)],[bounds.u_ub(1) bounds.u_ub(1)],'k--')
plot([t(1) t(end)],[bounds.u_lb(1) bounds.u_lb(1)],'k--')
title('$u_1$ - Force')
xlim([0 5])
subplot(4,2,6);hold on
plot(t,u(:,2),'k-')
plot([t(1) t(end)],[bounds.d_ub(1) bounds.d_ub(1)],'k--')
plot([t(1) t(end)],[bounds.d_lb(1) bounds.d_lb(1)],'k--')
title('$d_1$ - Road hieght change rate')
xlim([0 5])
subplot(4,2,7);hold on
plot(t,u(:,3),'k-')
plot(t,y(:,1),'r-')
plot(t,y(:,2),'b-')
legend('z_0','z_1','z_2','Orientation','horizontal')
title('$z_i$ - Heights')
xlim([0 5])
xlabel('Time (sec)')
subplot(4,2,8);hold on
plot(t,y(:,3),'k-')
title('$y_3$ - Sprung mass acceleration')
xlim([0 5])
xlabel('Time (sec)')
end