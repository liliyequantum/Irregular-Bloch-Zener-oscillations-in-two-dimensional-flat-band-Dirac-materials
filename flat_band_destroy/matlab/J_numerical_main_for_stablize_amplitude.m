clear
close all
maxNumCompThreads(6);

p_x_size=300;
p_y_size=500;
p_x=linspace(0,1/sqrt(3)*2*pi,p_x_size);
p_y=linspace(-1/3*2*pi,2/3*2*pi,p_y_size);
dp_x=p_x(2)-p_x(1);
dp_y=p_y(2)-p_y(1);

[p_x,p_y]=meshgrid(p_x,p_y);
p_x=p_x(:);
p_y=p_y(:);

global electric_field
global alpha_parameter
electric_field=0.0279;%4*pi/sqrt(3)/260;
alpha_parameter=0;

num_step=1e2;
t=linspace(0,36/electric_field,int64(36/electric_field*num_step));

% t=[0,0.01];
dt=t(2)-t(1);
t_max=t(end);

check_normalization = zeros(length(p_x),1);

x=zeros(length(p_x),6);
x(:,6)=1;

J_inter=zeros(length(t),1);
J_intra=zeros(length(t),1);

start = tic;
c=0;
for k = 1:1:length(t)-1

    T=[t(k),t(k+1)];
    tic

    %"""Runge Kutta fourth order"""
    x=RungeKutta4(T,x,p_x,p_y);%vector input, vector output

    %"""Current: J_p_inter and J_p_intra"""

    check_normalization= sum(x.*x,2);
    disp(strcat('normalization max: ',num2str(max(check_normalization)),' min:',num2str(min(check_normalization))))
         
    %compute J_matrix
    sin_2phi=2.*alpha_parameter./(1+(alpha_parameter).^2);
    cos_2phi=(1-(alpha_parameter).^2)./(1+(alpha_parameter).^2);
    f_k=-(1+2.*exp(-3./2.*1j.*p_y).*cos(sqrt(3)./2.*(p_x-electric_field.*t(k+1))));
    sin_theta_plus_p_y=sin(angle(f_k)+3./2.*p_y);
    cos_theta_plus_p_y=cos(angle(f_k)+3./2.*p_y);
    common_factor=-sqrt(3).*sin(sqrt(3)./2.*(p_x-electric_field.*t(k+1)));
    J_11 = common_factor.*cos_theta_plus_p_y;
    J_12 = common_factor.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y;
    J_13 = common_factor.*1j.*cos_2phi.*sin_theta_plus_p_y;
    J_23 = common_factor.*1j./sqrt(2).*sin_2phi.*sin_theta_plus_p_y;

    J_p_inter =J_p_inter_map(x,J_12,J_13,J_23) ;
    [J_p_intra,p_alpha_2_minus_beta_2] =J_p_intra_map_1(x,J_11);
    
    %integrate the first Brilloun hexagon region of momentum space
    J_inter(k+1,1)=sum(J_p_inter).*dp_x.*dp_y;
    J_intra(k+1,1)=sum(J_p_intra).*dp_x.*dp_y;
    
    if isnan(J_inter(k+1,1)) || isnan(J_intra(k+1,1))
        disp(['nan and k=',num2str(k)])
        save(['./data/Bloch_Oscillation_alpha_T3','_alpha_',num2str(alpha_parameter),...
    '_dp_x_',num2str(dp_x),'_dp_y_',num2str(dp_y),'_E_',num2str(electric_field),'_t_max_',num2str(t(k)),'_dt_',num2str(t(2)-t(1)),'.mat'], ...
    't','p_x','p_y','J_p_inter','J_p_intra','p_x_size','p_y_size','electric_field','J_inter','J_intra','alpha_parameter','x')
        stop
    end  
    disp(strcat("spend time: ", num2str(toc)," s")) 
    
    if k==int64(1/electric_field*num_step) || k==int64(2/electric_field*num_step) || k==int64(2.2/electric_field*num_step) || k==int64(2.4/electric_field*num_step) || k==int64(2.6/electric_field*num_step)...
      ||k==int64(3/electric_field*num_step) || k==int64(3.5/electric_field*num_step) || k==int64(4/electric_field*num_step) || k==int64(4.5/electric_field*num_step)...
      || k==int64(4.8/electric_field*num_step) || k==int64(5/electric_field*num_step) || k==int64(5.5/electric_field*num_step) || k==int64(6/electric_field*num_step)...
      || k==int64(6.5/electric_field*num_step) || k==int64(7/electric_field*num_step) || k==int64(7.2/electric_field*num_step) || k==int64(7.5/electric_field*num_step)...
      || k==int64(8/electric_field*num_step) || k==int64(8.5/electric_field*num_step) || k==int64(9/electric_field*num_step) || k==int64(9.5/electric_field*num_step) || k==int64(10/electric_field*num_step)...
      ||k==int64(10.5/electric_field*num_step) || k==int64(11/electric_field*num_step) || k==int64(12/electric_field*num_step) || k==int64(15/electric_field*num_step)...
      || k==int64(18/electric_field*num_step) || k==int64(20/electric_field*num_step) || k==int64(22/electric_field*num_step) || k==int64(24/electric_field*num_step)...
      || k==int64(28/electric_field*num_step) || k==int64(35/electric_field*num_step) || k==int64(40/electric_field*num_step) || k==int64(45/electric_field*num_step) || k==int64(50/electric_field*num_step)...
      ||k==int64(55/electric_field*num_step) || k==int64(60/electric_field*num_step) || k==int64(65/electric_field*num_step) || k==int64(70/electric_field*num_step)...
      || k==int64(75/electric_field*num_step) || k==int64(80/electric_field*num_step) || k==int64(85/electric_field*num_step) || k==int64(90/electric_field*num_step)...

        save(['./picture_data/J_t_max_',...
            num2str(t(k)),'_dt_',num2str(dt),'.mat'],'t','J_intra')
        
        J_p_inter = reshape(J_p_inter,[p_y_size,p_x_size]);
        J_p_intra = reshape(J_p_intra,[p_y_size,p_x_size]);
        p_alpha_2_minus_beta_2 = reshape(p_alpha_2_minus_beta_2,[p_y_size,p_x_size]);

        p_x=reshape(p_x,[p_y_size,p_x_size]);
        p_y=reshape(p_y,[p_y_size,p_x_size]);
        save(['./picture_data/p_alpha_2_minus_beta_2_t_max_',...
            num2str(t(k)),'_dt_',num2str(dt),'.mat'],'p_x','p_y','p_alpha_2_minus_beta_2')
        
        
        save(['./picture_data/intraband_t_max_',num2str(t(k)),...
            '_dt_',num2str(dt),'.mat'],'p_x','p_y','J_p_intra')
 
        
        save(['./data/Bloch_Oscillation_alpha_T3','_alpha_',num2str(alpha_parameter),...
    '_dp_x_',num2str(dp_x),'_dp_y_',num2str(dp_y),'_E_',num2str(electric_field),'_t_max_',num2str(t(k)),'_dt_',num2str(t(2)-t(1)),'.mat'], ...
    't','p_x','p_y','J_p_inter','J_p_intra','p_x_size','p_y_size','electric_field','J_inter','J_intra','alpha_parameter','x')
    end
    p_x=p_x(:);
    p_y=p_y(:);
    J_p_inter=J_p_inter(:);
    J_p_intra=J_p_intra(:);
end
disp(strcat("whole time: ", num2str(toc(start))," s")) 

save(['./data/Bloch_Oscillation_alpha_T3','_alpha_',num2str(alpha_parameter),...
    '_dp_x_',num2str(dp_x),'_dp_y_',num2str(dp_y),'_E_',num2str(electric_field),'_t_max_',num2str(t(end)),'_dt_',num2str(t(2)-t(1)),'.mat'], ...
    't','p_x','p_y','J_p_inter','J_p_intra','p_x_size','p_y_size','electric_field','J_inter','J_intra','alpha_parameter','x')



