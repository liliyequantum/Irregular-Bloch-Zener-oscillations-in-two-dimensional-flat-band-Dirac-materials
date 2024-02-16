clear
close all
color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

% k_x=1.5;
% k_y=2*pi/3+0.08;

k_x=0.5;
k_y=2*pi/3+0.07;
% 
% k_x=3;
% k_y=0.07;

% k_x=3;
% k_y=2*pi/3+0.07;

E=0.03;
t=-100:0.001:380;

factor=cos(sqrt(3)/2*(k_x-E*t));
energy_up=sqrt(1+4*factor.*(cos(3/2*k_y)+factor));
energy_down=-sqrt(1+4*factor.*(cos(3/2*k_y)+factor));

% [peak_down,loc_down]=findpeaks(energy_down);
% [peak_up,loc_up]=findpeaks(energy_up);

% t_peak_1=t(loc_down(1));
% t_peak_2=t(loc_down(2));
% t_2=((t_peak_1+t_peak_2)/2+t_peak_1)/2;
% t_1=t_peak_1-(t_2-t_peak_1);

% t_peak_1=t(loc_down(1));
% t_peak_2=t(loc_down(2));
% t_2=(((t_peak_1+t_peak_2)/2+t_peak_1)/2+t_peak_1)/2;
% t_1=t_peak_1-(t_2-t_peak_1);

% t_peak_1=t(loc_down(1));
% t_peak_2=t(loc_down(2));
% t_2=((((t_peak_1+t_peak_2)/2+t_peak_1)/2+t_peak_1)/2+t_peak_1)/2;
% t_1=t_peak_1-(t_2-t_peak_1);

% t_peak_1=t(loc_down(1));
% t_peak_2=t(loc_down(2));
% t_2=(((((t_peak_1+t_peak_2)/2+t_peak_1)/2+t_peak_1)/2+t_peak_1)/2+t_peak_1)/2;
% t_1=t_peak_1-(t_2-t_peak_1);


% factor_1=cos(sqrt(3)/2*(k_x-E*t_1));
% factor_2=cos(sqrt(3)/2*(k_x-E*t_2));
% energy_1=-sqrt(1+4*factor_1.*(cos(3/2*k_y)+factor_1));
% energy_2=sqrt(1+4*factor_2.*(cos(3/2*k_y)+factor_2));
% 
% S=(energy_2-energy_1)/(t_2-t_1)
% disp('gap')
% -max(peak_down)*2
% gap=-max(peak_down);
% S=2*E;
% probability=exp(-pi*gap^2/S)
% t_2=
plot(t,energy_up,'color',color_red,'LineWidth',2)
hold on
plot(t,energy_down,'color',color_blue,'LineWidth',2)
hold on
% plot([t_1,t_2],[energy_1,energy_2],'color',color_purple,'LineWidth',2)
% hold on
xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('Energy','Fontname', 'Times New Roman','FontSize',25)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlim([-100,380])
yticks([-3,0,3])
ylim([-3,3])
% load('E_0_12.mat')
% plot(t,alpha_2,'color',color_red,'LineWidth',2)
% hold on
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\alpha_{B}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

