clear
close all

% load('filter_hexagon_p_x_size_200_p_y_size_200.mat')%load filter_hexagon, p_x_size,p_y_size
% 
% % p_x=p_x+4*pi/sqrt(3);
% energy=sqrt(1+4.*cos(sqrt(3)/2.*p_x).*(cos(3/2.*p_y)+cos(sqrt(3)/2.*p_x)));
% a=4*pi/3/sqrt(3);
% b=0;
% energy_1=sqrt(1+4.*cos(sqrt(3)/2.*a).*(cos(3/2.*b)+cos(sqrt(3)/2.*a)))

% energy=energy.*filter_hexagon;
p_x_size=1200;
p_y_size=3000;
p_x=linspace(0,4*pi/3/sqrt(3),2*p_x_size);
p_y=linspace(-2*pi/3,2/3*2*pi,2.*p_y_size);
[p_x,p_y]=meshgrid(p_x,p_y);
p_x=p_x(:);
p_y=p_y(:);
energy_up=sqrt(1+4.*cos(sqrt(3)/2.*p_x).*(cos(3/2.*p_y)+cos(sqrt(3)/2.*p_x)));
energy_down=-sqrt(1+4.*cos(sqrt(3)/2.*p_x).*(cos(3/2.*p_y)+cos(sqrt(3)/2.*p_x)));

p_x=reshape(p_x,[2.*p_y_size,2.*p_x_size]);
p_y=reshape(p_y,[2.*p_y_size,2.*p_x_size]);

% energy_tra=reshape(energy,[2.*p_y_size,2.*p_x_size]);
energy_up=reshape(energy_up,[2.*p_y_size,2.*p_x_size]);
energy_down=reshape(energy_down,[2.*p_y_size,2.*p_x_size]);

figure(1)
surf(p_x,p_y,energy_up,'EdgeColor','none')
hold on
surf(p_x,p_y,energy_down,'EdgeColor','none')
hold on
colormap(parula(20));
colorbar;
xlabel('p_x','Fontname', 'Times New Roman','FontSize',16)
ylabel('p_y','Fontname', 'Times New Roman','FontSize',16)

figure(2)
plot(p_y(:,1200),energy_up(:,1200))
hold on
plot(p_y(:,1200),energy_down(:,1200))
hold on

figure(3)
surface(p_x,p_y,energy_up,'EdgeColor','none')
hold on
surface(p_x,p_y,energy_down,'EdgeColor','none')
hold on
colormap(parula(20));
colorbar;
xlabel('p_x','Fontname', 'Times New Roman','FontSize',16)
ylabel('p_y','Fontname', 'Times New Roman','FontSize',16)



% surf(p_x,p_y,energy)
% colormap(parula(20));

