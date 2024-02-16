clear
% close all
% color_blue = [31, 119, 180]/255;
% color_red = [214, 39, 40]/255;
% color_yellow = [255, 127, 14]/255;  
% color_purple = [148, 103, 189]/255;
% color_green = [44, 160, 44]/255;

dk_y=0:0.01:1;
% electric=[20,10,5,2.5,1.25,0.625,0.3125];
% k_x=1.5;
% t=-0.5:0.001:15;
% color=colorGradient([182, 33, 254]/255,[31, 209, 249]/255,length(electric));
% 
% gap=zeros(length(dk_y),length(electric));
% for j = 1:1:length(electric)
%     E=electric(j);
%     for i = 1:1:length(dk_y)
%         k_y=2*pi/3+dk_y(i);
% 
%         factor=cos(sqrt(3)/2*(k_x-E*t));
%         energy_down=-sqrt(1+4*factor.*(cos(3/2*k_y)+factor));
% 
%         [peak_down,loc_down]=findpeaks(energy_down);
% 
%         gap(i,j)=-max(peak_down)*2;
%     end
%     plot(dk_y,gap(:,j),'color',color(j,:),'LineWidth',2,'DisplayName',['E=',num2str(electric(j))])
%     hold on
% 
% end

analy_exp = 2*sqrt(1-cos(3/2*dk_y).^2);
plot(dk_y,analy_exp,'k-','LineWidth',2,'DisplayName','analytical expression')
hold on
xlabel('$\delta k_y$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
ylabel('gap','Fontname', 'Times New Roman','FontSize',16)
legend()
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)


