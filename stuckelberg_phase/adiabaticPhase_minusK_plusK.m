clear
close all

xpk = linspace(0,4*pi/3,1000);
xmk = linspace(0,2*pi/3,1000);

dky = 0.07;
E = 0.03;

ap = cos(1.5*dky);
am = -ap;
xipk = 2/sqrt(3)/E*sqrt(1+4*ap*cos(2*pi/3-xpk)+4*cos(2*pi/3-xpk).^2);
ximk = 2/sqrt(3)/E*sqrt(1+4*am*cos(pi/3-xmk)+4*cos(pi/3-xmk).^2);

color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;

figure()
plot(xpk/pi,xipk,'color',color_yellow,'LineWidth',3);hold on;
plot(xmk/pi,ximk,'color',color_purple,'LineWidth',3);hold on;
set(gcf,'Position',[100 100 1000 400])
text(1,90,'$+K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25,'color',color_yellow)
text(0.5,40,'$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25,'color',color_purple)
yticks([40,80,120])
xticks([0,1/3,2/3,3/3,4/3])
xticklabels({'0','1/3','2/3','1','4/3'})
ylabel('$f(x)$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
xlabel('$x(\pi)$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
axis tight