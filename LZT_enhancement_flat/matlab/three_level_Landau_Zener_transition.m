clear
close all

electric = 0.0336;
k_y = 0:0.001:0.4;
gap = 2*sqrt(1-cos(1.5*k_y).^2);
slope = 2*electric;
p = probability(gap,slope);

alpha_2 = p;
gamma_2 = 2*(1-sqrt(p)).*sqrt(p);
beta_2 = (1-sqrt(p)).^2;

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;


plot(k_y, alpha_2, 'color', color_red, 'LineWidth',2, 'DisplayName', '|\alpha|^2')
hold on
plot(k_y, gamma_2, 'color', color_green,'LineWidth',2, 'DisplayName', '|\gamma|^2')
hold on
plot(k_y, beta_2, 'color', color_blue,'LineWidth',2, 'DisplayName', '|\beta|^2')
hold on

legend('FontSize',16)
ylim([0,1])
xlabel('$k_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

save('./three_level_LZT_alpha_1_E_0_0336.mat','k_y','alpha_2','gamma_2','beta_2')
function p = probability(gap,slope)
    p = exp(-pi.*(gap./2).^2./slope);
end

