clear
close all

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

h=figure;
subplot(2,2,1)
load('./p_vs_k_y_alpha.mat')
surface(k_y, alpha, p,'EdgeColor','none')
xlabel('$\Delta\tilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
ylabel('$\alpha$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
xticks([0,0.1,0.2,0.3])
yticks([0,0.5,1.0])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',33)
colormap(parula(50));
cbh = colorbar;
cbh.Ticks = linspace(0, 1, 6) ; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(linspace(0, 1, 6)) ;  
cbh.FontSize = 30;
axis tight
hold on
plot3(0.1,1,100,'rp','MarkerSize',16, 'LineWidth',2)
hold on
% plot3(0.1,0,100,'r+','MarkerSize',16, 'LineWidth',2)
set(gcf,'Position',[10 10 1200 1000])
text(-0.1,1.1,'$(a)$','Interpreter','latex','FontSize',28)

subplot(2,2,2)
load('alpha_1_k_y_0_1_E_0_0336.mat')
plot(t,alpha_2_t,'color',color_red,'LineWidth',3,'Displayname','|\alpha_{k}|^{2}')
hold on
plot(t,gamma_2_t,'color',color_green,'LineWidth',3,'Displayname','|\gamma_{k}|^{2}')
hold on
plot(t,alpha_2_t + gamma_2_t,'color',color_purple,'LineWidth',3,'Displayname','|\alpha_{k}|^{2}+|\gamma_{k}|^{2}')
hold on
ylim([0,1])
yticks([0,0.25,0.5,0.75,1])
xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',35)
text(40,0.35,'$|\alpha_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_red)
text(40,0.6,'$|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_green)
text(40,0.85,'$|\alpha_{k}|^{2}+|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_purple)
set(gcf,'Position',[10 10 1200 1000])
text(-40,1.1,'$(b)$','Interpreter','latex','FontSize',28)

subplot(2,2,3)
% load('three_level_LZT_alpha_1_E_0_0336.mat')
k_y = 0:0.001:0.4;
E = 0.0336;
r = 1.5 * k_y.^2 / E;
p = exp( - pi  *r);
alpha_2 = p;
gamma_2 = 2 * (1 - sqrt(p)) .* sqrt(p);
plot(k_y, alpha_2, 'color', color_red, 'LineWidth',3, 'DisplayName', '|\alpha_{k}|^2')
hold on
plot(k_y, gamma_2, 'color', color_green,'LineWidth',3, 'DisplayName', '|\gamma_{k}|^2')
hold on
plot(k_y, alpha_2+gamma_2, 'color', color_purple,'LineWidth',3, 'DisplayName', '|\alpha_{k}|^2+|\gamma_{k}|^2')
hold on
plot([-0.05,0.4], [0.8,0.8], 'k--','LineWidth',1)
hold on
% legend('$|\alpha_{k}|^2$', '$|\gamma_{k}|^2$','$|\alpha_{k}|^2+|\gamma_{k}|^2$',...
%     'FontSize',25,'Interpreter','latex','location',[0.32,0.26,0.1,0.1],'box','off')
% text(0.04,0.85,'\rightarrow','FontSize',30,'Color','k','LineWidth',2)
text(0.2,0.4,'$|\alpha_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_red)
text(-0.04,0.4,'$|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_green)
text(0.12,0.7,'$|\alpha_{k}|^{2}+|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_purple)
ylim([0,1])
xlim([-0.05,k_y(end)])
xticks([0,0.2,0.4])
yticks([0,0.25,0.5,0.75,1])
xlabel('$\Delta\tilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',35)
set(gcf,'Position',[10 10 1200 1000])
text(-0.18,1.1,'$(c)$','Interpreter','latex','FontSize',28)

subplot(2,2,4)

load('./p_vs_alpha.mat')
plot(alpha, p_alpha, 'LineWidth',3,'color',color_red)
hold on
plot(alpha, p_gamma, 'LineWidth',3,'color',color_green)
hold on
% legend('$|\alpha_{k}|^2$', '$|\gamma_{k}|^2$',...
%     'FontSize',25,'Interpreter','latex','box','off','location','northwest')
text(0.25,0.3,'$|\alpha_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_red)
text(0.5,0.1,'$|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'color',color_green)
xlabel('$\alpha$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',28)
ylim([0,0.5])
xticks(0:0.5:1)
yticks([0,0.25,0.5])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',35)
set(gcf,'Position',[10 10 1200 1000])
text(-0.4,0.54,'$(d)$','Interpreter','latex','FontSize',28)

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')
saveas(gcf,'LZT_enhancement_flat.pdf')
