function plotSensitivity(app)

value = app.DisplaySensitivityDropDown.Value;

AtlasState = get_global_variable(app, 'AtlasState');
Sensitivity = get_global_variable(app, 'Sensitivity_Matrix');

% --- Create a new figure as a child of the app ---
fig = figure('Name', 'Atlas Viewer', ...
             'NumberTitle', 'off', ...
             'Color', [1 1 1], ...
             'Visible', 'on', ...
             'Units', 'normalized', ...
             'Position', [0.2 0.2 0.5 0.5], ...
             'HandleVisibility', 'callback', ...
             'CloseRequestFcn', @(src,evt) delete(src));
% So if app closes, this figure also closes
addlistener(app, 'ObjectBeingDestroyed', @(~,~) delete(fig));

% --- First axes (e.g., left plot)
ax1 = axes('Parent', fig, 'Position', [0.1 0.1 0.35 0.8]); 
title(ax1, 'Left Hemisphere');
hold(ax1, 'on');

% --- Second axes (e.g., right plot)
ax2 = axes('Parent', fig, 'Position', [0.55 0.1 0.35 0.8]);
title(ax2, 'Right Hemisphere');
hold(ax2, 'on');

% plot_intensity(ax, caxis_value,faces, brain_vertices, int_at_pos, rotation)

switch value
    case 'Head'
        faces       =   AtlasState.headsurf.mesh_reduced.faces;
        vertices    =   AtlasState.headsurf.mesh_reduced.vertices;
        Sens_value =    sum(Sensitivity.Adot_scalp(:,:,1),1);
    otherwise
        faces       =   AtlasState.pialsurf.mesh_reduced.faces;
        vertices    =   AtlasState.pialsurf.mesh_reduced.vertices;
        Sens_value = sum(Sensitivity.Adot(:,:,1),1);
end

plot_intensity(ax1, [-max(Sens_value), max(Sens_value)], faces,vertices,Sens_value,'L');
plot_intensity(ax2, [-max(Sens_value), max(Sens_value)], faces,vertices,Sens_value,'R');