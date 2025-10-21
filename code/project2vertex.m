function [HbO, HbR] = project2vertex(app,n_vertex_brain)
T = get_global_variable(app, 'T');
M = get_global_variable(app, 'M');
Conc = get_global_variable(app, 'Conc');

n_cond = size(T.T_HbO_brain,2);
tbasis = T.t_HbO_brain;
position = tbasis(:,:,1) == max(tbasis(:,:,1));

HbO = zeros(n_vertex_brain, n_cond);
HbR = zeros(n_vertex_brain, n_cond);

for i_cond = 1:n_cond
    HbO(M.mask_brain,i_cond) = Conc.intensity_HbO{i_cond}(:,position);
    HbR(M.mask_brain,i_cond) = Conc.intensity_HbR{i_cond}(:,position);
end