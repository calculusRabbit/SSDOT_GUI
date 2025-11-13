function [cfgOut, accepted] = userOptionGUI(cfgIn)
    % userOptionGUI - GUI for editing ssdot configuration parameters
    % Excludes 'device' as it's handled in the main GUI
    
    if nargin < 1 || isempty(cfgIn)
        cfgIn = struct( ...
            'threshold_brain', 5, ...
            'threshold_scalp', 5, ...
            'sigma_brain', 5, ...
            'sigma_scalp', 5, ...
            'spatially_regu', 0, ...
            'threshold_sensitivity', -2, ...
            'alpha', 1e-2, ...
            'beta', 0, ...
            'regularization', 1);
    end
    cfg = cfgIn; 
    accepted = false;

    fig = uifigure('Name', 'User Options', ...
                   'Position', [100 100 480 480], ...
                   'WindowStyle', 'modal', ...
                   'CloseRequestFcn', @(~,~) doClose());

    % Grid layout with 10 rows
    grid = uigridlayout(fig, [10 2]);
    grid.RowHeight   = {30, 30, 30, 30, 30, 30, 30, 30, 30, 44};
    grid.ColumnWidth = {200, '1x'};
    grid.Padding     = [12 12 12 12];
    grid.RowSpacing  = 8; 
    grid.ColumnSpacing = 10;

    % row 1: threshold_brain 
    uilabel(grid, 'Text', 'threshold_brain', 'HorizontalAlignment', 'right');
    efTB = uieditfield(grid, 'numeric', 'Value', cfg.threshold_brain, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 2: threshold_scalp
    uilabel(grid, 'Text', 'threshold_scalp', 'HorizontalAlignment', 'right');
    efTS = uieditfield(grid, 'numeric', 'Value', cfg.threshold_scalp, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 3: sigma_brain
    uilabel(grid, 'Text', 'sigma_brain', 'HorizontalAlignment', 'right');
    efSB = uieditfield(grid, 'numeric', 'Value', cfg.sigma_brain, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 4: sigma_scalp
    uilabel(grid, 'Text', 'sigma_scalp', 'HorizontalAlignment', 'right');
    efSS = uieditfield(grid, 'numeric', 'Value', cfg.sigma_scalp, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 5: treshold sensitivity
    uilabel(grid, 'Text', 'threshold_sensitivity', 'HorizontalAlignment', 'right');
    efTSens = uieditfield(grid, 'numeric', 'Value', cfg.threshold_sensitivity);

    % row 6: alpha
    uilabel(grid, 'Text', 'alpha', 'HorizontalAlignment', 'right');
    efAlpha = uieditfield(grid, 'numeric', 'Value', cfg.alpha, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 7: beta
    uilabel(grid, 'Text', 'beta', 'HorizontalAlignment', 'right');
    efBeta = uieditfield(grid, 'numeric', 'Value', cfg.beta, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 8: regularization
    uilabel(grid, 'Text', 'regularization', 'HorizontalAlignment', 'right');
    efReg = uieditfield(grid, 'numeric', 'Value', cfg.regularization, ...
        'LowerLimitInclusive', 'on', 'LowerLimit', 0);

    % row 9: spatially_regu checkbox ----
    uilabel(grid, 'Text', 'spatially_regu', 'HorizontalAlignment', 'right');
    cbSpatially = uicheckbox(grid, 'Text', 'Enable', ...
        'Value', logical(cfg.spatially_regu));

    % row 10: buttons
    btnCancel = uibutton(grid, 'Text', 'Cancel', 'ButtonPushedFcn', @(s,e) doClose());
    btnCancel.Layout.Row = 10;  
    btnCancel.Layout.Column = 1;

    btnOK = uibutton(grid, 'Text', 'OK', 'ButtonPushedFcn', @(s,e) onOK());
    btnOK.Layout.Row = 10;      
    btnOK.Layout.Column = 2;

    uiwait(fig);

    if accepted
        cfgOut = cfg;
    else
        cfgOut = cfgIn;
    end

    function onOK()
        % Read all values from GUI
        userInput.threshold_brain = efTB.Value;
        userInput.threshold_scalp = efTS.Value;
        userInput.sigma_brain = efSB.Value;
        userInput.sigma_scalp = efSS.Value;
        userInput.threshold_sensitivity = efTSens.Value;
        userInput.alpha = efAlpha.Value;
        userInput.beta = efBeta.Value;
        userInput.regularization = efReg.Value;
        userInput.spatially_regu = double(cbSpatially.Value);

        % Validate
        [ok, msg] = validateCfg(userInput);
        if ~ok
            uialert(fig, msg, 'Invalid Input'); 
            return;
        end

        cfg = userInput;
        accepted = true;
        doClose();
    end

    function doClose()
        if isvalid(fig)
            uiresume(fig);
            delete(fig);
        end
    end
end

function [ok, msg] = validateCfg(userInput)
    ok = true; 
    msg = "";

    % Non-negative fields
    nonnegFields = {'threshold_brain', 'threshold_scalp', ...
                    'sigma_brain', 'sigma_scalp', ...
                    'alpha', 'beta', 'regularization'};
    
    for k = 1:numel(nonnegFields)
        fld = nonnegFields{k};
        v = userInput.(fld);
        if ~isscalar(v) || ~isfinite(v) || v < 0
            ok = false; 
            msg = sprintf('%s must be a finite number >= 0.', fld); 
            return;
        end
    end

    % threshold_sensitivity: any finite scalar allowed (can be negative)
    v = userInput.threshold_sensitivity;
    if ~isscalar(v) || ~isfinite(v)
        ok = false; 
        msg = 'threshold_sensitivity must be a finite number.'; 
        return;
    end
    
    % spatially_regu this must be 0 or 1
    v = userInput.spatially_regu;
    if ~isscalar(v) || (v ~= 0 && v ~= 1)
        ok = false; 
        msg = 'spatially_regu must be 0 or 1.'; 
        return;
    end
end