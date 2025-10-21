function save_global_data(app, var_name, var_value)
% this func is to
% var should be a double value or a string

% checking: is it a string, or numbers



hFig =  app.UIFigure;
setappdata(hFig, var_name, var_value);

%check

