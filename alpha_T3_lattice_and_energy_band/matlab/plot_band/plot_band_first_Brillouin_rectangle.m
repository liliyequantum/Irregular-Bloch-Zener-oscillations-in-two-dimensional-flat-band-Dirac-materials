clear
close all

a=4*pi/3/sqrt(3);

p_x=linspace(-a,2.5*a,2000);
p_y=linspace(-sqrt(3)*a,1.5*sqrt(3)*a,1500);
dp_x=p_x(2)-p_x(1)
dp_y=p_y(2)-p_y(1)
[p_x,p_y]=meshgrid(p_x,p_y);
energy=sqrt(1+4.*cos(sqrt(3)/2.*p_x).*(cos(3/2.*p_y)+cos(sqrt(3)/2.*p_x)));

A=[4*pi/(3*sqrt(3)),2/3*2*pi];
B=[2*pi/(3*sqrt(3)),1/3*2*pi];
C=[4*pi/(3*sqrt(3)),0];
D=[2*pi/(3*sqrt(3)),-1/3*2*pi];

AA=[2*pi/sqrt(3),2/3*2*pi];
BB=[0,1/3*2*pi];
CC=[2*pi/sqrt(3),0];
DD=[0,-1/3*2*pi];

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

h=figure;

surface(p_x,p_y,energy,'EdgeColor','none')
hold on
plot3([AA(1),A(1),B(1),C(1),D(1),DD(1)],[AA(2),A(2),B(2),C(2),D(2),DD(2)],100*ones(1,6),'Color',color_red,'Linewidth', 3)
hold on
plot3([BB(1),B(1)],[BB(2),B(2)],100*ones(1,2),'Color',color_red,'Linewidth', 3)
hold on
plot3([CC(1),C(1)],[CC(2),C(2)],100*ones(1,2),'Color',color_red,'Linewidth', 3)
hold on
plot3([0,1.5*a,1.5*a,0,0],[sqrt(3)*a,sqrt(3)*a,-sqrt(3)*a/2,-sqrt(3)*a/2,sqrt(3)*a],99*ones(1,5),'k--','Linewidth', 2)
axis tight
% xlabel('p_x','Fontname', 'Times New Roman','FontSize',18)
% ylabel('p_y','Fontname', 'Times New Roman','FontSize',18)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',18)
width=480;
height=400;
set(gcf,'position',[10,10,width,height])
colormap(parula);
cbh = colorbar ; %Create Colorbar
% cbh.Ticks = [0,1,2,2.9] ; %Create 8 ticks from zero to 1
% cbh.TickLabels = {'0','1','2','2.9'} ;    %Replace the labels of these 8 ticks with the numbers 1 to 8

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')

saveas(gcf,'alpha_T3_band_rectangle.png')