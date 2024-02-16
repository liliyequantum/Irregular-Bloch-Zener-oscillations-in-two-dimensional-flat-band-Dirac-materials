clear
close all
maxNumCompThreads(6);

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

width_char = '0_5';
global electric_field
global alpha_parameter
electric_field=0.12;
alpha_parameter=0;

if length(width_char)==4
    num_step=5e2;
else
    num_step=1e2;
end

t=linspace(0,8/electric_field,int64(8/electric_field*num_step));
dt=t(2)-t(1);

c=0;
E_0=pi^2/4;
for k = 1:1:length(t)-1
%      if k==int64(1/electric_field*num_step) || k==int64(2/electric_field*num_step) || k==int64(2.2/electric_field*num_step) || k==int64(2.4/electric_field*num_step) || k==int64(2.6/electric_field*num_step)...
%       ||k==int64(3/electric_field*num_step) || k==int64(3.5/electric_field*num_step) || k==int64(4/electric_field*num_step) || k==int64(4.5/electric_field*num_step)...
%       || k==int64(4.8/electric_field*num_step) || k==int64(5/electric_field*num_step) || k==int64(5.5/electric_field*num_step) || k==int64(6/electric_field*num_step)...
%       || k==int64(6.5/electric_field*num_step) || k==int64(7/electric_field*num_step) || k==int64(7.2/electric_field*num_step) || k==int64(7.5/electric_field*num_step)...
%       || k==int64(8/electric_field*num_step) 
     if k==int64(7.5/electric_field*num_step)

        load(['./picture/E_0_12_small_k_area/picture_data/width_',width_char,...
        '/J_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.mat'])

        c=c+1;
        figure(c)
        plot(t(1:k),J_intra(1:k)./electric_field/E_0/3,'color',color_red,'LineWidth',6)
%         t = 0:0.01:9.5*pi;
%         y = sin(t);
%         plot(t,y,'color',color_red,'LineWidth',5)
        hold on
        set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
        axis off
        
%         saveas(gcf,['./picture/E_0_12_small_k_area/picture/width_',width_char,'/J_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.png'])
%         close
        
    end
end