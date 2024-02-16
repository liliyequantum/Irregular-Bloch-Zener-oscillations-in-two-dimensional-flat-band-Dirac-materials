clear
close all

color_yellow = [255, 127, 14]/255;  
color_purple = [148, 103, 189]/255;

k_0 = 4*pi/3/sqrt(3);
EE = [0.0725
0.06
0.0576
0.0554
0.0533
0.0514
0.0462
0.0448
0.0434
0.0386
0.0376
0.0366
0.0332
0.0317
0.0279
0.0072
0.0065
0.006
0.0059
0.0053
0.0614
0.0504
0.0487
0.0471
0.0344
0.0336
0.03
0.0415
0.0403
0.0392
0.0352
0.0294
0.0282
0.0261
0.0069
0.0068
0.0067
0.0062
0.0042
0.0035
0.003
];

stuckelberg_phase_vs_E = zeros(length(EE),1);
for j = 1:1:length(EE)
    
    E = EE(j);
    t = 0:0.01:2*k_0/E;
    dt = t(2)-t(1);
    k_x = k_0;

    K_y = 0.01:0.0001:sqrt(2/3*E); %r=1 to determin max(\Delta_\tilde{k}_y)
    phase_plus_k = zeros(length(K_y),1);
    for i=1:1:length(K_y)
        k_y = K_y(i);
        X_k_t = (k_x - E*t)*sqrt(3)/2;
        Y_k = k_y * 3/2;
        energy_p = sqrt(1+4*cos(X_k_t).*(cos(Y_k)+cos(X_k_t)));
        phase_plus_k(i) = mod(sum(energy_p)*dt,2*pi);
    end

    delta_k_y = K_y;
    delta = 3 * delta_k_y.^2 / 4 / E;
    phi_s = pi/4 + delta .* (log(delta) - 1) + angle(gamma(1 - 1j * delta));
    tilde_phi_s = mod(phi_s,2*pi);

    phase_stuckelberg = mod(phase_plus_k+tilde_phi_s',2*pi);
%     phase_stuckelberg = mod(phase_plus_k,2*pi);
    index = find(phase_stuckelberg>3*pi/2);
    phase_stuckelberg(index) = phase_stuckelberg(index)-2*pi;
   

    stuckelberg_phase_vs_E(j) = sum(phase_stuckelberg)/length(K_y);

end
figure(4)
scatter(1:1:20,stuckelberg_phase_vs_E(1:20),125,'o','LineWidth',2,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',color_yellow)
hold on
scatter(21:1:length(EE),stuckelberg_phase_vs_E(21:end),125,'v','LineWidth',2,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',color_yellow)
hold on
plot(1:1:length(EE),2*pi*ones(length(EE),1),'k-','linewidth',1)
hold on
plot(1:1:length(EE),3*pi/2*ones(length(EE),1),'k-','linewidth',1)
hold on
plot(1:1:length(EE),pi*ones(length(EE),1),'k-','linewidth',1)
hold on
plot(1:1:length(EE),pi/2*ones(length(EE),1),'k-','linewidth',1)
hold on
plot(1:1:length(EE),zeros(length(EE),1),'k-','linewidth',1)
hold on
plot([41,41],[-1,2*pi],'k-','linewidth',1)
set(gcf,'Position',[100 100 1000 400])
yticks([0,pi/2,pi,3*pi/2,2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xticks([1:4:length(EE)])
% axisHandle = gca;
% axisHandle.XAxis.TickLabelRotation = 90; % Rotate 90 counter clockwise.
% legend('$+K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'location','best')
ylim([-1,2*pi])
xlim([1,length(EE)])
ylabel('$\langle\Phi\rangle_{st}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
xlabel('$N$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
set(gcf,'Position',[100 100 1000 400])

% 
% figure(3)
% plot(K_y,mod(phase_plus_k+tilde_phi_s',2*pi),'color',color_yellow,'linewidth',3)
% hold on
% plot(K_y,2*pi*ones(length(K_y),1),'k-','linewidth',1)
% hold on
% plot(K_y,3*pi/2*ones(length(K_y),1),'k-','linewidth',1)
% hold on
% plot(K_y,pi*ones(length(K_y),1),'k-','linewidth',1)
% hold on
% plot(K_y,pi/2*ones(length(K_y),1),'k-','linewidth',1)
% hold on
% set(gcf,'Position',[100 100 1000 400])
% yticks([0,pi/2,pi,3*pi/2,2*pi])
% yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
% legend('$+K$','$-K$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30,'location','best')
% ylim([0,2*pi])
% xlim([K_y(1),K_y(end)])
% ylabel('$\Phi_{st}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
% xlabel('$\Delta\widetilde{k}_{y}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
% 

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
