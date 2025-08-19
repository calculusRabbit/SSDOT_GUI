function [HbO, HbR] = project2vertex(calculated_mat,n_vertex_brain)

n_cond = length(calculated_mat.Proc_data.stim);
position = calculated_mat.Proc_data.tbasis(:,:,1) == max(calculated_mat.Proc_data.tbasis(:,:,1));

HbO = zeros(n_vertex_brain, n_cond);
HbR = zeros(n_vertex_brain, n_cond);

for i_cond = 1:n_cond
    HbO(calculated_mat.M.mask_brain,i_cond) = calculated_mat.Conc.intensity_HbO{i_cond}(:,position);
    HbR(calculated_mat.M.mask_brain,i_cond) = calculated_mat.Conc.intensity_HbR{i_cond}(:,position);
end