clear
close all

electric=0.0336;
alpha = 0.3;

k_x=3;
k_y=0.1;

tspan=[0,100];
y0=[0,0,1];

t=tspan(1):0.01:tspan(2);
t=transpose(t);
y = RungeKutta4(t,y0,alpha,electric,k_x,k_y);

disp('t(2)-t(1)')
t(2)-t(1)

H_11=zeros(length(t),1);
H_12=zeros(length(t),1);
H_13=zeros(length(t),1);
H_21=zeros(length(t),1);
H_22=zeros(length(t),1);
H_23=zeros(length(t),1);
H_31=zeros(length(t),1);
H_32=zeros(length(t),1);
H_33=zeros(length(t),1);


for i = 1:1:length(t)
    H=Hamiltonian(t(i),alpha,electric,k_x,k_y);
    H_11(i)=H(1,1);
    H_12(i)=H(1,2);
    H_13(i)=H(1,3);
    H_21(i)=H(2,1);
    H_22(i)=H(2,2);
    H_23(i)=H(2,3);
    H_31(i)=H(3,1);
    H_32(i)=H(3,2);
    H_33(i)=H(3,3);
end

error_1=diff(y(:,1))./diff(t)-(H_11(2:end).*y(2:end,1)+H_12(2:end).*y(2:end,2)+H_13(2:end).*y(2:end,3));
error_2=diff(y(:,2))./diff(t)-(H_21(2:end).*y(2:end,1)+H_22(2:end).*y(2:end,2)+H_23(2:end).*y(2:end,3));
error_3=diff(y(:,3))./diff(t)-(H_31(2:end).*y(2:end,1)+H_32(2:end).*y(2:end,2)+H_33(2:end).*y(2:end,3));

disp('max(real(error_1)) and max(imag(error_1))')
max(abs(real(error_1)))
max(abs(imag(error_1)))
disp('max(real(error_2)) and max(imag(error_2))')
max(abs(real(error_2)))
max(abs(imag(error_2)))
disp('max(real(error_3)) and max(imag(error_3))')
max(abs(real(error_3)))
max(abs(imag(error_3)))

alpha_2_t=abs(y(:,1)).^2;
gamma_2_t=abs(y(:,2)).^2;
beta_2_t=abs(y(:,3)).^2;

check_norm=alpha_2_t+gamma_2_t+beta_2_t;
disp('normalization max and min')
max(check_norm)
min(check_norm)

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

% figure(1)
% plot(t,alpha_2_t,'color',color_red,'LineWidth',2)
% hold on
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% ylim([0,1])
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\alpha_{k}(t)|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% 
% figure(2)
% plot(t,gamma_2_t,'color',color_green,'LineWidth',2)
% hold on
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% ylim([0,1])
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\gamma_{k}(t)|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

% figure(3)
% plot(t,beta_2_t,'color',color_blue,'LineWidth',2)
% hold on
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% % ylim([0,1])
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\beta_{k}(t)|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

% figure(4)
% plot(t,2*alpha_2_t + gamma_2_t -1,'color',color_purple,'LineWidth',2)
% hold on
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% % ylim([0,1])
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\alpha_{k}(t)|^{2}-|\beta_{k}(t)|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

% figure(5)
% plot(t,alpha_2_t + gamma_2_t,'color',color_purple,'LineWidth',2)
% hold on
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% ylim([0,1])
% xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% ylabel('$|\alpha_{k}|^{2}+|\gamma_{k}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

figure(6)
plot(t,alpha_2_t,'color',color_red,'LineWidth',2,'Displayname','|\alpha_{k}(t)|^{2}')
hold on
plot(t,gamma_2_t,'color',color_green,'LineWidth',2,'Displayname','|\gamma_{k}(t)|^{2}')
hold on
plot(t,alpha_2_t + gamma_2_t,'color',color_purple,'LineWidth',2,'Displayname','|\alpha_{k}(t)|^{2}+|\gamma_{k}(t)|^{2}')
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
ylim([0,1])
xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
text(60,0.35,'|\alpha_{k}(t)|^{2}|','Fontname', 'Times New Roman','FontSize',16,'color',color_red)
text(60,0.6,'|\gamma_{k}(t)|^{2}|','Fontname', 'Times New Roman','FontSize',16,'color',color_green)
text(60,0.85,'|\alpha_{k}(t)|^{2}+|\gamma_{k}(t)|^{2}','Fontname', 'Times New Roman','FontSize',16,'color',color_purple)

save(['alpha_',num2str(alpha),'_k_y_0_1_E_0_0336.mat'],'t','alpha_2_t','gamma_2_t')
function energy=band(t,electric,k_x,k_y)
factor=cos(sqrt(3)/2.*(k_x-electric.*t));
energy=sqrt(1+4.*factor.*(cos(3/2.*k_y)+factor));
end

function H=Hamiltonian(t,alpha,electric,k_x,k_y)
sin_2phi=2.*alpha./(1+(alpha).^2);
cos_2phi=(1-(alpha).^2)./(1+(alpha).^2);
ep=band(t,electric,k_x,k_y);
C_0=1/sqrt(2)*sqrt(3)*electric*sin(3/2*(-k_y))*sin(sqrt(3)/2*(k_x-electric*t))/ep^2;
H=-1j*[ep+1/sqrt(2)*C_0*cos_2phi,    C_0*sin_2phi,                 1/sqrt(2)*C_0*cos_2phi;
       C_0*sin_2phi,                 -sqrt(2)*C_0*cos_2phi,        C_0*sin_2phi;
    1/sqrt(2)*C_0*cos_2phi,           C_0*sin_2phi,                -ep+1/sqrt(2)*C_0*cos_2phi];
end

function y = RungeKutta4(t,y0,alpha,electric,k_x,k_y)
   
    y=zeros(length(t),3);
    y(1,:)=y0;
    step_size=t(2)-t(1);
    
    for i = 1:1:length(t)-1
        k_1=f(t(i),y(i,:),alpha,electric,k_x,k_y);
        k_2=f(t(i)+0.5.*step_size,y(i,:)+0.5.*step_size.*k_1,alpha,electric,k_x,k_y);
        k_3=f(t(i)+0.5.*step_size,y(i,:)+0.5.*step_size.*k_2,alpha,electric,k_x,k_y);
        k_4=f(t(i)+step_size,y(i,:)+step_size.*k_3,alpha,electric,k_x,k_y);

        y(i+1,:)=y(i,:)+1/6.*step_size.*(k_1+2.*k_2+2.*k_3+k_4);
    end
        
end

function y = f(t,y0,alpha,electric,k_x,k_y)
    
   y = zeros(1,3);
   H=Hamiltonian(t,alpha,electric,k_x,k_y);
   y(1) = H(1,1).*y0(1) + H(1,2).*y0(2) + H(1,3).*y0(3);
   y(2) = H(2,1).*y0(1) + H(2,2).*y0(2) + H(2,3).*y0(3);
   y(3) = H(3,1).*y0(1) + H(3,2).*y0(2) + H(3,3).*y0(3);
 
end


