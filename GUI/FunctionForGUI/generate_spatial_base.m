function generate_spatial_base(app)
%% Calculate G
atlasViewer = get_global_variable(app, 'AtlasState');
Sensitivity_Matrix = get_global_variable(app, 'Sensitivity_Matrix');

% Vu we need to get threshold_brain and _scalp from user input
cfg.threshold_brain = 5;
cfg.threshold_scalp = 20;
cfg.sigma_brain = 5;
cfg.sigma_scalp = 20;
% Vu we should save all the user input data in cfg
hFig =  app.UIFigure;
setappdata(hFig, 'cfg', cfg);

G = Make_G_matrix(atlasViewer, Sensitivity_Matrix.M, ...
    cfg.threshold_brain, cfg.threshold_scalp, cfg.sigma_brain, cfg.sigma_scalp, ...
    Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);

hFig =  app.UIFigure;
setappdata(hFig, 'G', G);

disp('calculate G finished');