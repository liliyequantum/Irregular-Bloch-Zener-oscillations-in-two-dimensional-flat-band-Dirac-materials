clear
close all

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

% electric_field=0:0.001:15;
electric_field=0.01:0.001:0.2;
gap=[2,1.8,1.5,1.3,1,0.8,2^(-1),2^(-2),2^(-3),2^(-4)];
% gap=[2^(-1),2^(-2),2^(-3),2^(-4)];

color=colorGradient([182, 33, 254]/255,[31, 209, 249]/255,length(gap));
probability=zeros(length(gap),length(electric_field));

for i =1:1:length(gap)
probability(i,:)=exp(-pi*(gap(i)/2).^2./(2*electric_field));
plot(electric_field,probability(i,:),'LineWidth',2,'color',color(i,:))
hold on
end

% plot([5,5],[0,1],'k--','linewidth',1)
% hold on
% plot([0,15],[0.732,0.732],'k--','linewidth',1)
% hold on
% 
% plot([1.26,1.26],[0,1],'k--','linewidth',1)
% hold on
% plot([0,15],[0.3,0.3],'k-','linewidth',1)
% hold on

plot([0.03,0.03],[0,1],'k--','linewidth',1)
hold on
plot([0.12,0.12],[0,1],'k--','linewidth',1)
hold on

legend('gap=2','gap=1.8','gap=1.5','gap=1.3','gap=1','gap=0.8','gap=0.5',...
    'gap=0.25','gap=0.125','gap=0.0625','location','best',...
    'Fontname','Times New Roman','FontSize',14)

ylim([0,1])
xlim([electric_field(1),electric_field(end)])
yticks([0.3,0.6,0.8,1.0])
xticks([0.03,0.06,0.09,0.12,0.15,0.18])

xlabel('$\tilde{E}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
