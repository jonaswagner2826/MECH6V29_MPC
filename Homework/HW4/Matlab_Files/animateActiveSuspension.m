function animateActiveSuspension(t,u,y,road_x,road_z,pos,scale)
z1 = y(:,1); % unsprung mass height
z2 = y(:,2); % sprung mass height
zmf = 3;     % exaggerate response for better visualization
speedUp = 10;% Speed up animation by showing every x frames
f1 = figure;
for i=1:speedUp:length(t)
    plotsusp([u(i,3), z1(i)*zmf, z2(i)*zmf, t(i)],road_x,road_z,pos(i),scale,f1);
%     refresh
end
end