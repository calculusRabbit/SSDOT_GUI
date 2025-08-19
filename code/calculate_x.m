function [abs_x,Ax_y] = calculate_x(calculated_mat)
n_brain = size(calculated_mat.H_brain, 2);
index = 1:n_brain;

b_bs = calculated_mat.b(index);
b = calculated_mat.b;
abs_x = b_bs*b_bs';
Y_prime = [calculated_mat.H_brain, calculated_mat.H_scalp, calculated_mat.OD_SS, calculated_mat.OD_drift]*b';
Y_prime_Y = Y_prime - calculated_mat.Y;
Ax_y = Y_prime_Y'*Y_prime_Y;
end