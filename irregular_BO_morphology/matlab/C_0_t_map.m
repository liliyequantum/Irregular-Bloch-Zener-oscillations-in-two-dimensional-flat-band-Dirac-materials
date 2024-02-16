function  C_0_t = C_0_t_map(t,p_x,p_y)
global electric_field
numerator=sqrt(3).*electric_field.*sin(3/2.*(-p_y)).*sin(sqrt(3)./2.*(p_x-electric_field.*t));
epsilon_p_t = epsilon_p_t_map(t,p_x,p_y);
denominator= epsilon_p_t.^2;
C_0_t=1/sqrt(2).*numerator./denominator;
end