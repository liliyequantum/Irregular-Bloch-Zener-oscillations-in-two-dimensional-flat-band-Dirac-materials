clear
close all

color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;

k_0 = 4*pi/3/sqrt(3);
E = 0.031;

t = 0:0.001:2*k_0/E;
dt = t(2)-t(1);
k_x = k_0;
K_y = 0.01:0.0001:sqrt(2/3*E);

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

delta_k_y = K_y;
delta = 3 * delta_k_y.^2 / 4 / E;
phi_s = pi/4 + delta .* (log(delta) - 1) + angle(gamma(1 - 1j * delta));
tilde_phi_s = mod(phi_s,2*pi);

phase_stuckelberg = mod(phase_plus_k+tilde_phi_s',2*pi);
stuckelberg_phase_vs_E = sum(phase_stuckelberg)/length(K_y)
    
figure(1)
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
yticks([-pi/2,0,pi,2*pi])
yticklabels({'-\pi/2','0','\pi','2\pi'})
legend('$+K$','$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20,'location','best')
ylim([-pi/2,2*pi])
xlim([K_y(1),K_y(end)])
ylabel('$\zeta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
xlabel('$\Delta\widetilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)

figure(2)
plot(delta_k_y,tilde_phi_s,'linewidth',3)
hold on
plot(delta_k_y,2*pi*ones(length(delta_k_y),1),'k-','linewidth',1)
hold on
plot(delta_k_y,pi*ones(length(delta_k_y),1),'k-','linewidth',1)
hold on

yticks([0,pi,2*pi])
yticklabels({'0','\pi','2\pi'})
ylim([0,2*pi])
xlim([delta_k_y(1),delta_k_y(end)])
ylabel('$\widetilde{\varphi}_{s}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
xlabel('$\Delta\widetilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)

figure(3)
phi_st_plus_k = mod(phase_plus_k+tilde_phi_s',2*pi);
index = find(phi_st_plus_k>3*pi/2);
phi_st_plus_k(index) = phi_st_plus_k(index)-2*pi;

plot(K_y,phi_st_plus_k,'color',color_yellow,'linewidth',3)
hold on
plot(K_y,mod(phase_minus_k+tilde_phi_s',2*pi),'color',color_purple,'linewidth',3)
hold on
plot(K_y,2*pi*ones(length(K_y),1),'k-','linewidth',1)
hold on
plot(K_y,3*pi/2*ones(length(K_y),1),'k-','linewidth',1)
hold on
plot(K_y,pi*ones(length(K_y),1),'k-','linewidth',1)
hold on
plot(K_y,pi/2*ones(length(K_y),1),'k-','linewidth',1)
hold on
set(gcf,'Position',[100 100 1000 400])
text(0.08,0.5,'$+K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25,'color',color_yellow)
text(0.033,pi+0.5,'$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25,'color',color_purple)
yticks([-pi/2,0,pi/2,pi,3*pi/2,2*pi])
yticklabels({'-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'})
% legend('$+K$','$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'location','best')
ylim([-pi/2,2*pi])
xlim([K_y(1),K_y(end)])
ylabel('$\Phi_{st}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
xlabel('$\Delta\widetilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)

function [f] = gamma(z)
% GAMMA  Gamma function valid in the entire complex plane.
%        Accuracy is 15 significant digits along the real axis
%        and 13 significant digits elsewhere.
%        This routine uses a superb Lanczos series
%        approximation for the complex Gamma function.
%
%        z may be complex and of any size.
%        Also  n! = prod(1:n) = gamma(n+1)
%
%usage: [f] = gamma(z)
%       
%tested on versions 6.0 and 5.3.1 under Sun Solaris 5.5.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931-944
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%            W. J. Cody "An Overview of Software Development for Special
%            Functions", 1975
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%
%Paul Godfrey
%pgodfrey@intersil.com
%http://winnie.fit.edu/~gabdo/gamma.txt
%Sept 11, 2001

siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
   z(p)=-z(p);
end

% 15 sig. digits for 0<=real(z)<=171
% coeffs should sum to about g*g/2+23/24

g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];

%Num Recipes used g=5 with 7 terms
%for a less effective approximation

z=z-1;
zh =z+0.5;
zgh=zh+g;
%trick for avoiding FP overflow above z=141
zp=zgh.^(zh*0.5);

ss=0.0;
for pp=size(c,1)-1:-1:1
    ss=ss+c(pp+1)./(z+pp);
end

%sqrt(2Pi)
sq2pi=  2.5066282746310005024157652848110;
f=(sq2pi*(c(1)+ss)).*((zp.*exp(-zgh)).*zp);

f(z==0 | z==1) = 1.0;

%adjust for negative real parts
if ~isempty(p)
   f(p)=-pi./(zz(p).*f(p).*sin(pi*zz(p)));
end

%adjust for negative poles
p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
   f(p)=Inf;
end

f=reshape(f,siz);

% return
end
