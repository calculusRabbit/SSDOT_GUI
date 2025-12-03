function generate_spatial_base(app, ssdot)
%% Calculate G
% atlasViewer = get_global_variable(app, 'AtlasState');
atlasViewer = ssdot.getVar('AtlasState');

% Sensitivity_Matrix = get_global_variable(app, 'Sensitivity_Matrix');
Sensitivity_Matrix = ssdot.getVar('Sensitivity_Matrix');

ssdot.setVar('cfg', cfg);

G = Make_G_matrix(atlasViewer, Sensitivity_Matrix.M, ...
    cfg.threshold_brain, cfg.threshold_scalp, cfg.sigma_brain, cfg.sigma_scalp, ...
    Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);

% hFig =  app.UIFigure;
% setappdata(hFig, 'G', G);
ssdot.setVar('G', G);

disp('calculate G finished');