clear
close all
% maxNumCompThreads(6);

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

width_char = '0_5';
global electric_field
global alpha_parameter
electric_field=0.12;
alpha_parameter=0;

if length(width_char)==4
    num_step=5e2;
else
    num_step=1e3;
end

numPeriod =  20; % Bloch period
t=linspace(0,(8*numPeriod)/electric_field,int64((8*numPeriod)/electric_field*num_step));
tmp = 1:0.5:(8*numPeriod);
selectTime = int64(tmp/electric_field*num_step);
dt=t(2)-t(1);

c=318;
E_0=pi^2/4;


for k = 1329167%1:1:length(t)-1
     if k==selectTime(c)
            
        load(['./picture/E_0_12_small_k_area/picture_data/width_',width_char,...
            '/p_alpha_2_minus_beta_2_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.mat'])
        
        
        figure(c)
        surface(p_x_1,p_y_1,p_alpha_2_minus_beta_2_1,'EdgeColor','none')
        hold on
        yticks([-0.4,-0.2,0,0.2,0.4])

        colormap(parula(20));
        colorbar;

        set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)

        axis tight
        saveas(gcf,['./picture/E_0_12_small_k_area/picture/width_',width_char,'/p_alpha_2_minus_beta_2_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.png'])
        close

        c=c+1;
    end
end