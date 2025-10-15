function loadSens(app)

folderPath = uigetdir(pwd, 'Select AtlasViewer Folder which contains fw folder');

if isequal(folderPath, 0)
    disp('User canceled folder selection.');
    return;
end

disp(['Selected folder: ', folderPath]);

mlActAuto = []; % we need to do something here later
atlasViewer = get_global_variable(app, 'AtlasState');
Sensitivity_Matrix = Get_A_dot(folderPath,7,mlActAuto,1,atlasViewer);

% calculated_mat.Adot = Sensitivity_Matrix.Adot;
% calculated_mat.Adot_scalp = Sensitivity_Matrix.Adot_scalp;
% calculated_mat.E = Sensitivity_Matrix.E;
% calculated_mat.channels = Sensitivity_Matrix.channels;
% calculated_mat.shortSepChLst = Sensitivity_Matrix.shortSepChLst;
% if cfg.spatially_regu
%     calculated_mat.ll = Sensitivity_Matrix.ll;
%     calculated_mat.ll_scalp = Sensitivity_Matrix.ll_scalp;
% end

hFig =  app.UIFigure;
setappdata(hFig, 'Sensitivity_Matrix', Sensitivity_Matrix);

disp('Sensitivity_Matrix loaded and stored in GUI.');