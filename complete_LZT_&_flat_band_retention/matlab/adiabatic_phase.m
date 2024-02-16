clear
close all

k_0 = 4*pi/3/sqrt(3);
E = 0.0504;

t = 0:0.01:2*k_0/E;
dt = t(2)-t(1);
k_x = k_0;
K_y = 0.01:0.001:0.15;

phase_plus_k = zeros(length(K_y),1);
phase_minus_k = zeros(length(K_y),1);
for i=1:1:length(K_y)
    k_y = K_y(i);
    X_k_t = (k_x - E*t)*sqrt(3)/2;
    Y_k = k_y * 3/2;
    energy_p = sqrt(1+4*cos(X_k_t).*(cos(Y_k)+cos(X_k_t)));
    phase_plus_k(i) = mod(sum(energy_p)*dt,2*pi);
    
    k_y = K_y(i) + 2*pi/3;
    X_k_t = (k_x/2 - E*t/2)*sqrt(3)/2;
    Y_k = k_y * 3/2;
    energy_p = sqrt(1+4*cos(X_k_t).*(cos(Y_k)+cos(X_k_t)));
    phase_minus_k(i) = mod(sum(energy_p)*dt,2*pi);

end
plot(K_y,phase_plus_k,'linewidth',3)
hold on
plot(K_y,phase_minus_k,'linewidth',3)
hold on
plot(K_y,2*pi*ones(length(K_y),1),'k-','linewidth',1)
hold on
plot(K_y,pi*ones(length(K_y),1),'k-','linewidth',1)
hold on
% text(-0.001,6.6,'2\pi','FontSize',20)
% text(-0.000000001,3.3,'\pi','FontSize',20)
yticks([0,pi,2*pi])
yticklabels({'0','\pi','2\pi'})
legend('$+K$','$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20,'location','best')
ylim([0,2*pi])
xlim([K_y(1),K_y(end)])
ylabel('$\zeta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
xlabel('$\Delta\widetilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
