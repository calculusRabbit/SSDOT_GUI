function T = Make_T_matrix(dA,tbasis)

T_HbO_brain = squeeze(dA(:,:,1));
T_HbR_brain = squeeze(dA(:,:,2)); 
t_HbO_brain = squeeze(tbasis(:,:,1));
t_HbR_brain = squeeze(tbasis(:,:,2)); 


T.T_HbO_brain = T_HbO_brain;
T.T_HbR_brain = T_HbR_brain;

T.t_HbO_brain = t_HbO_brain;
T.t_HbR_brain = t_HbR_brain;

T.T_HbO_scalp = T_HbO_brain;
T.T_HbR_scalp = T_HbR_brain; 