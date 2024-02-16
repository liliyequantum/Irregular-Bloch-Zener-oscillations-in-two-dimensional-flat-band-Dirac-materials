function plot_J_p_intra

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

A=[4*pi/(3*sqrt(3)),2/3*2*pi];
B=[2*pi/(3*sqrt(3)),1/3*2*pi];
C=[4*pi/(3*sqrt(3)),0];
D=[2*pi/(3*sqrt(3)),-1/3*2*pi];

AA=[2*pi/sqrt(3),2/3*2*pi];
BB=[0,1/3*2*pi];
CC=[2*pi/sqrt(3),0];
DD=[0,-1/3*2*pi];

global electric_field
global alpha_parameter
electric_field=0.0279;%4*pi/sqrt(3)/260;
alpha_parameter=0;

num_step=1e2;
t=linspace(0,36/electric_field,int64(36/electric_field*num_step));
dt=t(2)-t(1);

c=0;
for k = 1:1:length(t)-1

if k==int64(1/electric_field*num_step) || k==int64(2/electric_field*num_step) || k==int64(2.2/electric_field*num_step) || k==int64(2.4/electric_field*num_step) || k==int64(2.6/electric_field*num_step)...
  ||k==int64(3/electric_field*num_step) || k==int64(3.5/electric_field*num_step) || k==int64(4/electric_field*num_step) || k==int64(4.5/electric_field*num_step)...
  || k==int64(4.8/electric_field*num_step) || k==int64(5/electric_field*num_step) || k==int64(5.5/electric_field*num_step) || k==int64(6/electric_field*num_step)...
  || k==int64(6.5/electric_field*num_step) || k==int64(7/electric_field*num_step) || k==int64(7.2/electric_field*num_step) || k==int64(7.5/electric_field*num_step)...
  || k==int64(8/electric_field*num_step) || k==int64(8.5/electric_field*num_step) || k==int64(9/electric_field*num_step) || k==int64(9.5/electric_field*num_step) || k==int64(10/electric_field*num_step)...
  ||k==int64(10.5/electric_field*num_step) || k==int64(11/electric_field*num_step) || k==int64(12/electric_field*num_step) || k==int64(15/electric_field*num_step)...
  || k==int64(18/electric_field*num_step) || k==int64(20/electric_field*num_step) || k==int64(22/electric_field*num_step) || k==int64(24/electric_field*num_step)...
  || k==int64(28/electric_field*num_step) || k==int64(35/electric_field*num_step) || k==int64(40/electric_field*num_step) || k==int64(45/electric_field*num_step) || k==int64(50/electric_field*num_step)...
  ||k==int64(55/electric_field*num_step) || k==int64(60/electric_field*num_step) || k==int64(65/electric_field*num_step) || k==int64(70/electric_field*num_step)...
  || k==int64(75/electric_field*num_step) || k==int64(80/electric_field*num_step) || k==int64(85/electric_field*num_step) || k==int64(90/electric_field*num_step)...

load(['./picture_data/intraband_t_max_',num2str(t(k)),...
            '_dt_',num2str(dt),'.mat'])
c=c+1;
figure(c)
surface(p_x,p_y,J_p_intra,'EdgeColor','none')
hold on
plot3([AA(1),A(1),B(1),C(1),D(1),DD(1)],[AA(2),A(2),B(2),C(2),D(2),DD(2)],100*ones(1,6),'Color',color_red,'Linewidth', 1)
hold on
plot3([BB(1),B(1)],[BB(2),B(2)],100*ones(1,2),'Color',color_red,'Linewidth', 1)
hold on
plot3([CC(1),C(1)],[CC(2),C(2)],100*ones(1,2),'Color',color_red,'Linewidth', 1)
hold on
%         title('intra band')
colormap(parula(20));
colorbar;
% xlabel('k_x','Fontname', 'Times New Roman','FontSize',16)
% ylabel('k_y','Fontname', 'Times New Roman','FontSize',16)
xticks([0.5,1.5,2.5,3.5])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
axis tight
saveas(gcf,['./picture/intraband_t_max_',...
    num2str(t(k)),'_dt_',num2str(dt),'.png'])
close
end
end