function plot_Hb(app, faces, brain_vertices, HbO, HbR)
caxis_value = [-1 1]*1e-9; % Vu need user input

for i_cond = 1:size(HbO,2)
    % --- Create a new figure as a child of the app ---
    fig = figure('Name', 'HbO', ...
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

    % HbO
    intensity_HbO = HbO(:,i_cond);
    plot_intensity(ax1, caxis_value, faces, brain_vertices, intensity_HbO, 'L')
    plot_intensity(ax2, caxis_value, faces, brain_vertices, intensity_HbO, 'R')
    % HbR
    % --- Create a new figure as a child of the app ---
    fig2 = figure('Name', 'HbR', ...
                 'NumberTitle', 'off', ...
                 'Color', [1 1 1], ...
                 'Visible', 'on', ...
                 'Units', 'normalized', ...
                 'Position', [0.2 0.2 0.5 0.5], ...
                 'HandleVisibility', 'callback', ...
                 'CloseRequestFcn', @(src,evt) delete(src));
    % So if app closes, this figure also closes
    addlistener(app, 'ObjectBeingDestroyed', @(~,~) delete(fig2));
    
    % --- First axes (e.g., left plot)
    ax1 = axes('Parent', fig2, 'Position', [0.1 0.1 0.35 0.8]); 
    title(ax1, 'Left Hemisphere');
    hold(ax1, 'on');
    % --- Second axes (e.g., right plot)
    ax2 = axes('Parent', fig2, 'Position', [0.55 0.1 0.35 0.8]);
    title(ax2, 'Right Hemisphere');
    hold(ax2, 'on');


    intensity_HbR = HbR(:,i_cond);
    plot_intensity(ax1, caxis_value, faces,brain_vertices, intensity_HbR, 'L')
    plot_intensity(ax2, caxis_value, faces,brain_vertices, intensity_HbR, 'R')
end
end