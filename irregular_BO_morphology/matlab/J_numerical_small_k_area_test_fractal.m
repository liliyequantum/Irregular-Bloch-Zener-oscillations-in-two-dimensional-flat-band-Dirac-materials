clear
close all
% maxNumCompThreads(10);

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

A=[4*pi/(3*sqrt(3)),2/3*2*pi];
B=[2*pi/(3*sqrt(3)),1/3*2*pi];
C=[4*pi/(3*sqrt(3)),0];
D=[2*pi/(3*sqrt(3)),-1/3*2*pi];

AA=[2*pi/sqrt(3),2/3*2*pi];
BB=[0,1/3*2*pi];
CC=[2*pi/sqrt(3),0];
DD=[0,-1/3*2*pi];

p_x_size=500;
p_y_size=500;
% p_x_1=linspace(1,2.4,p_x_size);
p_x_1=linspace(2.3,2.8,p_x_size);
% p_x_2=linspace(0,1/sqrt(3)*2*pi,p_x_size);

width_y = 0.5;
width_char = '0_5';
% p_y_1=linspace(2*pi/3-width_y,2*pi/3+width_y,p_y_size);
p_y_1=linspace(-width_y,+width_y,p_y_size);
% p_y_2=linspace(-width_y,+width_y,p_y_size);

dp_x=p_x_1(2)-p_x_1(1)
dp_y=p_y_1(2)-p_y_1(1)

[p_x_1,p_y_1]=meshgrid(p_x_1,p_y_1);
% [p_x_2,p_y_2]=meshgrid(p_x_2,p_y_2);

p_x_1=p_x_1(:);
p_y_1=p_y_1(:);
% p_x_2=p_x_2(:);
% p_y_2=p_y_2(:);

global electric_field
global alpha_parameter
electric_field=0.12;
alpha_parameter=0;

% time step number
if length(width_char)==4
    num_step=5e2;
else
    num_step=1e3;
end    

numPeriod = 20; % Bloch period
t=linspace(0,(8*numPeriod)/electric_field,int64((8*numPeriod)/electric_field*num_step));
tmp = 1:0.5:(8*numPeriod);
selectTime = int64(tmp/electric_field*num_step);

% t=[0,0.01];
dt=t(2)-t(1)
t_max=t(end);

check_normalization_1 = zeros(length(p_x_1),1);
% check_normalization_2 = zeros(length(p_x_2),1);

x_1=zeros(length(p_x_1),6);
x_1(:,6)=1;
% x_2=zeros(length(p_x_2),6);
% x_2(:,6)=1;

J_inter=zeros(length(t),1);
J_intra=zeros(length(t),1);

% start = tic;
c=1;
checkNormalizationMax = zeros(length(t)-1,1);
checkNormalizationMin = zeros(length(t)-1,1);
for k = 1:1:length(t)-1

    T=[t(k),t(k+1)];
%     tic

    %"""Runge Kutta fourth order"""
    x_1=RungeKutta4(T,x_1,p_x_1,p_y_1);%vector input, vector output
%     x_2=RungeKutta4(T,x_2,p_x_2,p_y_2);%vector input, vector output
    
    %"""Current: J_p_inter and J_p_intra"""
    checkNormalizationMax(k, 1)= max(sum(x_1.*x_1,2));
%     checkNormalizationMax(k, 2)= max(sum(x_2.*x_2,2));
    checkNormalizationMin(k, 1)= min(sum(x_1.*x_1,2));
%     checkNormalizationMin(k, 2)= min(sum(x_2.*x_2,2));
     
    %compute J_matrix_1
    sin_2phi=2.*alpha_parameter./(1+(alpha_parameter).^2);
    cos_2phi=(1-(alpha_parameter).^2)./(1+(alpha_parameter).^2);
    f_k_1=-(1+2.*exp(-3./2.*1j.*p_y_1).*cos(sqrt(3)./2.*(p_x_1-electric_field.*t(k+1))));
    sin_theta_plus_p_y_1=sin(angle(f_k_1)+3./2.*p_y_1);
    cos_theta_plus_p_y_1=cos(angle(f_k_1)+3./2.*p_y_1);
    common_factor_1=-sqrt(3).*sin(sqrt(3)./2.*(p_x_1-electric_field.*t(k+1)));
    J_11_1 = common_factor_1.*cos_theta_plus_p_y_1;
    J_12_1 = common_factor_1.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y_1;
    J_13_1 = common_factor_1.*1j.*cos_2phi.*sin_theta_plus_p_y_1;
    J_23_1 = common_factor_1.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y_1;

    J_p_inter_1 =J_p_inter_map(x_1,J_12_1,J_13_1,J_23_1) ;
    [J_p_intra_1,p_alpha_2_minus_beta_2_1] =J_p_intra_map_1(x_1,J_11_1);
    
    %J_maxtrix_2
    %compute J_matrix
      
%     f_k_2=-(1+2.*exp(-3./2.*1j.*p_y_2).*cos(sqrt(3)./2.*(p_x_2-electric_field.*t(k+1))));
%     sin_theta_plus_p_y_2=sin(angle(f_k_2)+3./2.*p_y_2);
%     cos_theta_plus_p_y_2=cos(angle(f_k_2)+3./2.*p_y_2);
%     common_factor_2=-sqrt(3).*sin(sqrt(3)./2.*(p_x_2-electric_field.*t(k+1)));
%     J_11_2 = common_factor_2.*cos_theta_plus_p_y_2;
%     J_12_2 = common_factor_2.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y_2;
%     J_13_2 = common_factor_2.*1j.*cos_2phi.*sin_theta_plus_p_y_2;
%     J_23_2 = common_factor_2.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y_2;
% 
%     J_p_inter_2 =J_p_inter_map(x_2,J_12_2,J_13_2,J_23_2) ;
%     [J_p_intra_2,p_alpha_2_minus_beta_2_2] =J_p_intra_map_1(x_2,J_11_2);
    
    %integrate the first Brilloun hexagon region of momentum space
    J_inter(k+1,1)=sum(J_p_inter_1).*dp_x.*dp_y;
    J_intra(k+1,1)=sum(J_p_intra_1).*dp_x.*dp_y;
    
    
    if isnan(J_inter(k+1,1)) || isnan(J_intra(k+1,1))
        disp(['nan and k=',num2str(k)])
        save(['./data/Bloch_Oscillation_alpha_T3','_alpha_',num2str(alpha_parameter),...
    '_dp_x_',num2str(dp_x),'_dp_y_',num2str(dp_y),'_E_',num2str(electric_field),'_t_max_',num2str(t(k)),'_dt_',num2str(t(2)-t(1)),'.mat'], ...
    't','p_x_1','p_y_1','J_p_inter_1','J_p_intra_1','p_x_size','p_y_size','electric_field','J_inter','J_intra','alpha_parameter','x_1')
        stop
    end  
%     disp(strcat("spend time: ", num2str(toc)," s")) 
    
    if k==selectTime(c) 
        
        save(['./picture/E_0_12_small_k_area/picture_data/width_',width_char,...
            '/J_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.mat'],'t','J_intra')

        
        J_p_inter_1 = reshape(J_p_inter_1,[p_y_size,p_x_size]);
        J_p_intra_1 = reshape(J_p_intra_1,[p_y_size,p_x_size]);
        p_alpha_2_minus_beta_2_1 = reshape(p_alpha_2_minus_beta_2_1,[p_y_size,p_x_size]);
        
%         J_p_inter_2 = reshape(J_p_inter_2,[p_y_size,p_x_size]);
%         J_p_intra_2 = reshape(J_p_intra_2,[p_y_size,p_x_size]);
%         p_alpha_2_minus_beta_2_2 = reshape(p_alpha_2_minus_beta_2_2,[p_y_size,p_x_size]);
        
        p_x_1=reshape(p_x_1,[p_y_size,p_x_size]);
        p_y_1=reshape(p_y_1,[p_y_size,p_x_size]);
        
%         p_x_2=reshape(p_x_2,[p_y_size,p_x_size]);
%         p_y_2=reshape(p_y_2,[p_y_size,p_x_size]);

        save(['./picture/E_0_12_small_k_area/picture_data/width_',width_char,...
            '/p_alpha_2_minus_beta_2_t_max_',num2str(t(k)),'_dt_',num2str(dt),'.mat'],...
            'p_x_1','p_y_1','p_alpha_2_minus_beta_2_1')

        disp(k/(length(t)-1))
        c = c+1;
    end
    p_x_1=p_x_1(:);
    p_y_1=p_y_1(:);
%     p_x_2=p_x_2(:);
%     p_y_2=p_y_2(:);
    J_p_inter_1=J_p_inter_1(:);
    J_p_intra_1=J_p_intra_1(:);
%     J_p_inter_2=J_p_inter_2(:);
%     J_p_intra_2=J_p_intra_2(:);
end
% disp(strcat("whole time: ", num2str(toc(start))," s")) 
disp(['check normalization max: ',num2str(max(max(checkNormalizationMax)))])
disp(['check normalization min: ',num2str(min(min(checkNormalizationMin)))])

