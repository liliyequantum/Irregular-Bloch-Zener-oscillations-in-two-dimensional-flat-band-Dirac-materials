clear
close all

electric=0.0336;
alpha = 0:0.01:1;

k_x=3;
k_y=0.1;

tspan=[0,100];
y0=[0,0,1];

p_alpha = zeros(length(alpha),1);
p_gamma = zeros(length(alpha),1);
p_beta = zeros(length(alpha),1);
for j = 1:1:length(alpha)
    
    tic
    t=tspan(1):0.01:tspan(2);
    t=transpose(t);
    y = RungeKutta4(t,y0,alpha(j),electric,k_x,k_y);
    alpha_2_t=abs(y(:,1)).^2;
    gamma_2_t=abs(y(:,2)).^2;
    beta_2_t=abs(y(:,3)).^2;
    p_alpha(j) = alpha_2_t(end);
    p_gamma(j) = gamma_2_t(end);
    p_beta(j) = beta_2_t(end);
    toc
    
end

close
color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;
% load('./p_vs_alpha.mat')
plot(alpha, p_alpha, 'LineWidth',3,'color',color_red)
hold on
plot(alpha, p_gamma, 'LineWidth',3,'color',color_green)
hold on
plot(alpha, p_beta, 'LineWidth',3,'color',color_blue)
hold on
plot(alpha, p_alpha+p_gamma+p_beta, 'LineWidth',3,'color',color_purple)
hold on

legend('$|\alpha_{k}|^2$', '$|\gamma_{k}|^2$', '$|\beta_{k}|^2$','total',...
    'FontSize',20,'Interpreter','latex','box','off','location','northwest')

xlabel('$\alpha$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
ylabel('$P$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
ylim([0,1])
xticks(0:0.2:1)
yticks([0,0.2,0.4,0.6,0.8,1.0])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
% axis tight

save('./p_vs_alpha.mat','alpha','p_alpha','p_gamma')

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


