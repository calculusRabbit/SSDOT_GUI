function generate_spatial_base(app)

%% Prompt the user for a mask threshold value
prompt = {'Enter threshold of sensitivity'};
dlgtitle = 'Input Threshold';
dims = [1 35];
definput = {'-2'};  % Default value

answer = inputdlg(prompt, dlgtitle, dims, definput);

% If user presses Cancel, exit
if isempty(answer)
    return;
end

% Convert input from string to numeric
mask_threshold = str2double(answer{1});

% Validate the input
if isnan(mask_threshold)
    uialert(app.UIFigure, 'Invalid input. Please enter a numeric value.', 'Error');
    return;
end
%% make a mask
atlasViewer = get_global_variable(app, 'AtlasState');
M = Make_mask(mask_threshold, Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);





brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

mask_brain = M.mask_brain;
mask_scalp = M.mask_scalp;

brain_vertices_masked = brain_vertices(mask_brain,:);
scalp_vertices_masked = scalp_vertices(mask_scalp,:);

[brain_vertices_new] = down_sample_vertices(brain_vertices_masked, threshold_brain);
[scalp_vertices_new] = down_sample_vertices(scalp_vertices_masked, threshold_scalp);


G_brain = make_kernel_matrix_gpu(brain_vertices_new, brain_vertices_masked,sigma_brain);
G_scalp = make_kernel_matrix_gpu(scalp_vertices_new, scalp_vertices_masked,sigma_scalp);