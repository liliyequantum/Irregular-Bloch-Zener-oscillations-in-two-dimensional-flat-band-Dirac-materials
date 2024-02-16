clear
close all

electric=0.0336;
alpha = 0:0.01:1;

k_x=3;
k_y=0.01:0.01:0.3;
% k_y=[-0.3:0.01:-0.01,0.01:0.01:0.3];

tspan=[0,100];
y0=[0,0,1];

p = zeros(length(alpha),length(k_y));
for j = 1:1:length(alpha)
    for i = 1:1:length(k_y)
        tic
        t=tspan(1):0.01:tspan(2);
        t=transpose(t);
        y = RungeKutta4(t,y0,alpha(j),electric,k_x,k_y(i));
        beta_2_t=abs(y(:,3)).^2;
        temp = 1-beta_2_t;
        p(j,i) = temp(end);
        toc
    end
end

% p = p(:,1:30);
[k_y,alpha] = meshgrid(k_y,alpha);

surface(k_y, alpha, p,'EdgeColor','none')
xlabel('k_y','Fontname', 'Times New Roman','FontSize',16)
ylabel('\alpha','Fontname', 'Times New Roman','FontSize',16)
xticks([0,0.1,0.2,0.3])
yticks([0,0.2,0.4,0.6,0.8,1.0])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',16)
colormap(parula(50));
cbh = colorbar;
cbh.Ticks = linspace(0, 1, 6) ; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(linspace(0, 1, 6)) ;  
axis tight

save('./p_vs_k_y_alpha.mat','k_y','alpha','p')
% saveas(gcf,['./picture/intraband_t_max_',...
%     num2str(t(k)),'_dt_',num2str(dt),'.png'])

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


