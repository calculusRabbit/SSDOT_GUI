function var = get_global_variable(app, var_name)

% Get the AtlasState structure stored earlier
hFig = app.UIFigure;
var = getappdata(hFig, var_name);

% Check if it exists
if isempty(var_name)
    uialert(app.UIFigure, ['Please load ', var_name,' data first.'], 'Data Missing');
    return;
end