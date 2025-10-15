function loadAV(app, event)

% load atlasviewer state file
[file, path] = uigetfile;
if isequal(file,0)
    disp('User canceled file selection');
    return;
end

AtlasState = load(fullfile(path,file));

% we need to save this AtlasState for global use; Ask Vu to
% also save data somewhere; we need OD data for each channel;
% we also need time basis
hFig =  app.UIFigure;
setappdata(hFig, 'AtlasState', AtlasState);

disp('AtlasState loaded and stored in GUI.');