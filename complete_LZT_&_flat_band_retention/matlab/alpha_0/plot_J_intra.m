function plot_J_intra

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
electric_field=0.03;
alpha_parameter=0;

num_step=1e3;
t=linspace(0,36/electric_field,int64(36/electric_field*num_step));
dt=t(2)-t(1);

c=0;
for k = 1:1:length(t)-1

% if  k==int64(2/electric_field*num_step) || k==int64(2.2/electric_field*num_step) || k==int64(2.4/electric_field*num_step) || k==int64(2.6/electric_field*num_step)...
%   ||k==int64(3/electric_field*num_step) || k==int64(3.5/electric_field*num_step) || k==int64(4/electric_field*num_step) || k==int64(4.5/electric_field*num_step)...
%   || k==int64(4.8/electric_field*num_step) || k==int64(5/electric_field*num_step) || k==int64(5.5/electric_field*num_step) || k==int64(6/electric_field*num_step)...
%   || k==int64(6.5/electric_field*num_step) || k==int64(7/electric_field*num_step) || k==int64(7.2/electric_field*num_step) || k==int64(7.5/electric_field*num_step)...
%   || k==int64(8/electric_field*num_step) || k==int64(8.5/electric_field*num_step) || k==int64(9/electric_field*num_step) || k==int64(9.5/electric_field*num_step) || k==int64(10/electric_field*num_step)...
%   ||k==int64(10.5/electric_field*num_step) || k==int64(11/electric_field*num_step) || k==int64(12/electric_field*num_step) || k==int64(15/electric_field*num_step)...
%   || k==int64(18/electric_field*num_step) || k==int64(20/electric_field*num_step) || k==int64(22/electric_field*num_step) || k==int64(24/electric_field*num_step)...
%   || k==int64(28/electric_field*num_step) || k==int64(35/electric_field*num_step) || k==int64(40/electric_field*num_step) || k==int64(45/electric_field*num_step) || k==int64(50/electric_field*num_step)...
%   ||k==int64(55/electric_field*num_step) || k==int64(60/electric_field*num_step) || k==int64(65/electric_field*num_step) || k==int64(70/electric_field*num_step)...
%   || k==int64(75/electric_field*num_step) || k==int64(80/electric_field*num_step) || k==int64(85/electric_field*num_step) || k==int64(90/electric_field*num_step)...
if k==int64(35/electric_field*num_step)
t(k)

load(['./picture_data/J_t_max_',...
            num2str(t(k)),'_dt_',num2str(dt),'.mat'])
E_0=pi^2/4;
c=c+1;
figure(c)
%         plot(t(1:k),J_inter(1:k)./electric_field/E_0/3,'color',color_blue,'LineWidth',2)
%         hold on
plot(t(1:k),J_intra(1:k)./electric_field/E_0/3,'color',color_red,'LineWidth',4)
hold on
%         plot(t(1:k),(J_intra(1:k)+J_inter(1:k))./electric_field/E_0/3,'color',color_yellow,'LineWidth',2)
%         hold on

if exist('picture','dir')==0
    mkdir('picture');
end
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% ylabel('$\tilde{J}/\tilde{E}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
xlim([t(1),t(k)])
ylim([-6,9])
xticks([0,500,1000])
axis off

set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
saveas(gcf,['./picture/J_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.png'])
% close

end

end