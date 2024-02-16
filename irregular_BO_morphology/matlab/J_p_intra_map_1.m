% function [J_p_intra,J_p_upper,J_p_flat] = J_p_intra_map(x,t,p_x,p_y)
function [J_p_intra,p_alpha_2_minus_beta_2] = J_p_intra_map_1(x,J_11)

alpha_2=x(:,4).^2+x(:,1).^2;
% gamma_2=x(:,5).^2+x(:,2).^2;
beta_2=x(:,6).^2+x(:,3).^2;
% J_p_intra=J_11.*(2.*alpha_2+gamma_2);
p_alpha_2_minus_beta_2=alpha_2-beta_2;
J_p_intra=J_11.*(p_alpha_2_minus_beta_2+1);

% J_p_upper=-(p_x-electric_field.*t)./epsilon_p_t.*(2.*alpha_2);
% J_p_flat=-(p_x-electric_field.*t)./epsilon_p_t.*gamma_2;
end
