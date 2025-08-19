function CNR = calculate_cnr(calculated_mat, cfg, Conc)

cfg.with_scalp = 0;
position = calculated_mat.Proc_data.tbasis(:,:,1) == max(calculated_mat.Proc_data.tbasis(:,:,1));
intensity = Conc.intensity_HbO{1};
intensity = intensity(:,position);

%% image contrast :  average of image contrast over FVHM
FVHM_index = find(intensity > 0.5*max(intensity));
CNR.Contrast = max(intensity(FVHM_index)) - min(intensity(FVHM_index));
%% image CNR : contrast/average of noise over FVHM
CNR.image_CNR = CNR.Contrast / std(intensity);
%% image localization error: 
At_file = 'atlasViewer.mat';
atlasViewer = load([cfg.av_data_path, At_file]);
brain_vertices = atlasViewer.fwmodel.mesh.vertices; 
error = 0;
sum_x = 0;
axes_order = [2,1,3];

    
for i = 1:length(FVHM_index)
    index = FVHM_index(i);
    x = intensity(index)-min(intensity(FVHM_index));
    index_brain = calculated_mat.M.mask_brain(index);
    r = sqrt((brain_vertices(index_brain,axes_order(1)) - cfg.center(1)).^2 + ...
        (brain_vertices(index_brain,axes_order(2)) - cfg.center(2)).^2 + ...
        (brain_vertices(index_brain,axes_order(3)) - cfg.center(3)).^2);
    error = error + x*r;
    sum_x = sum_x + x;
end
if length(FVHM_index) == 1
    CNR.error = r;
else
    CNR.error = error/sum_x;
end
end
