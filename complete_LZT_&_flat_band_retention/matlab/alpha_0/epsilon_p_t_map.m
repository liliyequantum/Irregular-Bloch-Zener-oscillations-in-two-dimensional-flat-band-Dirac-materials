function epsilon_p_t = epsilon_p_t_map(t,p_x,p_y)
global electric_field
epsilon_p_t=sqrt(1+4.*cos(sqrt(3)/2.*(p_x-electric_field.*t)).*(cos(3/2.*p_y)+cos(sqrt(3)/2.*(p_x-electric_field.*t))));
index=find(abs(epsilon_p_t)<=1e-5);
epsilon_p_t(index)=1e-5;
end
