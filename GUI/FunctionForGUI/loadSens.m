function loadSens(app, ssdot)

folderPath = uigetdir(pwd, 'Select AtlasViewer Folder which contains fw folder');

if isequal(folderPath, 0)
    disp('User canceled folder selection.');
    return;
end

disp(['Selected folder: ', folderPath]);

mlActAuto = []; % we need to do something here later
% atlasViewer = get_global_variable(app, 'AtlasState');
atlasViewer = ssdot.getVar('AtlasState');
% wavelength = [760 850]; % Vu we need users to input this or we get this information from snirf
run1 = ssdot.getVar('RunList');
run1 = run1(1);
wavelength = run1.GetWls();
%% Prompt the user for a mask threshold value
% prompt = {'Enter threshold of sensitivity'};
% dlgtitle = 'Input Threshold';
% dims = [1 35];
% definput = {'-2'};  % Default value

% answer = inputdlg(prompt, dlgtitle, dims, definput);
cfg = ssdot.getVar('cfg');
mask_threshold = cfg.threshold_sensitivity;

% If user presses Cancel, exit
% if isempty(answer)
%     return;
% end



% Convert input from string to numeric
%mask_threshold = str2double(answer{1});

% Validate the input
if isnan(mask_threshold)
    uialert(app.UIFigure, 'Invalid input. Please enter a numeric value.', 'Error');
    return;
end


% wavelength = [760 850]; % Vu please also get users to input this or we get this information from snirf
spatially_regu = 0; % Vu: get user input, either 0-no or 1-yes
Sensitivity_Matrix = Get_A_dot(folderPath,7,mlActAuto,spatially_regu,atlasViewer, wavelength, mask_threshold);

% hFig =  app.UIFigure;
% setappdata(hFig, 'Sensitivity_Matrix', Sensitivity_Matrix);
ssdot.setVar('Sensitivity_Matrix', Sensitivity_Matrix);

M = Make_mask(mask_threshold, Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);
% save_global_data(app,'M',M);
ssdot.setVar('M', M);

A = Make_A_matrix(Sensitivity_Matrix.Adot, Sensitivity_Matrix.Adot_scalp, Sensitivity_Matrix.E,  M);
% save_global_data(app,'A',A);
ssdot.setVar('A', A);

disp('calculate A finished');