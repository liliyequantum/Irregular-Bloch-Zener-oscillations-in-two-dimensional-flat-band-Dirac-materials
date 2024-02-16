clear
close all

% electric=4*pi/sqrt(3)/0.005
electric=0.0642;

% k_x=1.5;
% k_y=2*pi/3+0.08;

% k_x=0.5;
% k_y=2*pi/3+0.07;
% p=exp(-pi*(0.87/2)^2/electric)
% p=exp(-pi*(0.299/2)^2/2/electric)

k_x=3;
k_y=0.03;

% k_x=3;
% k_y=2*pi/3+0.15;

tspan=[0,1000];
% y0=[-0.0187788398103013-0.823907379601095i,-0.565350686732724-0.0346782067061930i]; %t=10
% y0 = [0.783333862616546 - 0.251956832174371i,-0.193682592623508 + 0.534221732711940i]; %t=15
% y0 = [-0.782188353041936 - 0.262113458247905i,-0.119076778736774 - 0.552538349911596i]; %t=18
y0=[0,1];

% opts = odeset('RelTol',1e-5,'AbsTol',1e-8);
% [t,y] = ode45(@(t,y) odefcn(t,y,electric,k_x,k_y), tspan, y0,opts);

t=tspan(1):0.01:tspan(2);
t=transpose(t);
y = RungeKutta4(t,y0,electric,k_x,k_y);

disp('t(2)-t(1)')
t(2)-t(1)

H_11=zeros(length(t),1);
H_12=zeros(length(t),1);
H_21=zeros(length(t),1);
H_22=zeros(length(t),1);

for i = 1:1:length(t)
    H=Hamiltonian(t(i),electric,k_x,k_y);
    H_11(i)=H(1,1);
    H_12(i)=H(1,2);
    H_21(i)=H(2,1);
    H_22(i)=H(2,2);
end

error_1=diff(y(:,1))./diff(t)-(H_11(2:end).*y(2:end,1)+H_12(2:end).*y(2:end,2));
error_2=diff(y(:,2))./diff(t)-(H_21(2:end).*y(2:end,1)+H_22(2:end).*y(2:end,2));
disp('max(real(error_1)) and max(imag(error_1))')
max(abs(real(error_1)))
max(abs(imag(error_1)))
disp('max(real(error_2)) and max(imag(error_2))')
max(abs(real(error_2)))
max(abs(imag(error_2)))

alpha_2_t=abs(y(:,1)).^2;
beta_2_t=abs(y(:,2)).^2;

check_norm=alpha_2_t+beta_2_t;
disp('normalization max and min')
max(check_norm)
min(check_norm)

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;


figure(1)
plot(t,alpha_2_t,'color',color_red,'LineWidth',2)
hold on
% plot(t,beta_2_t,'color',color_blue,'LineWidth',2)
% hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% legend('|\alpha|^{2}','|\beta|^{2}')
% legend('|\alpha|^{2}')
ylim([0,1])
xlabel('$\tilde{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
ylabel('$|\alpha_{B}|^{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
% t=15:0.01:30;
% energy=band(t,electric,k_x,k_y);
% figure(2)
% plot(t,energy,'color',color_red,'LineWidth',2)
% hold on
% plot(t,-energy,'color',color_red,'LineWidth',2)
% hold on

set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)

% function dydt = odefcn(t,y,electric,k_x,k_y)
%   dydt = zeros(2,1);
%   H=Hamiltonian(t,electric,k_x,k_y);
%   dydt(1) = H(1,1).*y(1) + H(1,2).*y(2);
%   dydt(2) = H(2,1).*y(1) + H(2,2).*y(2);
% end

function energy=band(t,electric,k_x,k_y)
factor=cos(sqrt(3)/2.*(k_x-electric.*t));
energy=sqrt(1+4.*factor.*(cos(3/2.*k_y)+factor));
end

function H=Hamiltonian(t,electric,k_x,k_y)

ep=band(t,electric,k_x,k_y);
C_0=1/sqrt(2)*sqrt(3)*electric*sin(3/2*(-k_y))*sin(sqrt(3)/2*(k_x-electric*t))/ep^2;
H=-1j*[ep+1/sqrt(2)*C_0,1/sqrt(2)*C_0;
    1/sqrt(2)*C_0,-ep+1/sqrt(2)*C_0];
end

function y = RungeKutta4(t,y0,electric,k_x,k_y)
   
    y=zeros(length(t),2);
    y(1,:)=y0;
    step_size=t(2)-t(1);
    
    for i = 1:1:length(t)-1
        k_1=f(t(i),y(i,:),electric,k_x,k_y);
        k_2=f(t(i)+0.5.*step_size,y(i,:)+0.5.*step_size.*k_1,electric,k_x,k_y);
        k_3=f(t(i)+0.5.*step_size,y(i,:)+0.5.*step_size.*k_2,electric,k_x,k_y);
        k_4=f(t(i)+step_size,y(i,:)+step_size.*k_3,electric,k_x,k_y);

        y(i+1,:)=y(i,:)+1/6.*step_size.*(k_1+2.*k_2+2.*k_3+k_4);
    end
        
end

function y = f(t,y0,electric,k_x,k_y)
    
   y = zeros(1,2);
   H=Hamiltonian(t,electric,k_x,k_y);
   y(1) = H(1,1).*y0(1) + H(1,2).*y0(2);
   y(2) = H(2,1).*y0(1) + H(2,2).*y0(2);
 
end


