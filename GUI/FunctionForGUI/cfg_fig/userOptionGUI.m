function [cfgOut, accepted] = userOptionGUI(cfgIn)


    if nargin<1 || isempty(cfgIn)
        cfgIn = struct( ...
            'threshold_brain', 5, ...
            'threshold_scalp', 5, ...
            'sigma_brain', 5, ...
            'sigma_scalp', 5, ...
            'threshold_sensitivity', -2);
    end
    cfg = cfgIn; 
    accepted = false;


    fig = uifigure('Name','UserOptionGUI', ...
                   'Position',[100 100 460 250], ...  
                   'WindowStyle','modal', ...
                   'CloseRequestFcn', @(~,~) doClose());

    
    g = uigridlayout(fig,[6 2]);
    g.RowHeight   = {30,30,30,30,30,44};
    g.ColumnWidth = {200,'1x'};
    g.Padding     = [12 12 12 12];
    g.RowSpacing  = 8; 
    g.ColumnSpacing = 10;

    % ---- labels + fields ----
    uilabel(g,'Text','threshold_brain','HorizontalAlignment','right');
    efTB = uieditfield(g,'numeric','Value',cfg.threshold_brain, ...
        'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','threshold_scalp','HorizontalAlignment','right');
    efTS = uieditfield(g,'numeric','Value',cfg.threshold_scalp, ...
        'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','sigma_brain','HorizontalAlignment','right');
    efSB = uieditfield(g,'numeric','Value',cfg.sigma_brain, ...
        'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','sigma_scalp','HorizontalAlignment','right');
    efSS = uieditfield(g,'numeric','Value',cfg.sigma_scalp, ...
        'LowerLimitInclusive','on','LowerLimit',0);

    
    uilabel(g,'Text','threshold_sensitivity','HorizontalAlignment','right');
    efTSens = uieditfield(g,'numeric','Value',cfg.threshold_sensitivity);

    % ---- buttons ----
    btnCancel = uibutton(g,'Text','Cancel','ButtonPushedFcn',@(s,e) doClose());
    btnCancel.Layout.Row = 6;  btnCancel.Layout.Column = 1;

    btnOK = uibutton(g,'Text','OK','ButtonPushedFcn',@(s,e) onOK());
    btnOK.Layout.Row = 6;      btnOK.Layout.Column = 2;

    uiwait(fig);

    if accepted, cfgOut = cfg; else, cfgOut = cfgIn; end

   
    function onOK()
        % Read values
        P.threshold_brain = efTB.Value;
        P.threshold_scalp = efTS.Value;
        P.sigma_brain = efSB.Value;
        P.sigma_scalp = efSS.Value;
        P.threshold_sensitivity = efTSens.Value;

        % Validate
        [ok,msg] = validateCfg(P);
        if ~ok
            uialert(fig,msg,'Invalid input'); 
            return;
        end

        cfg = P;
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

function [ok,msg] = validateCfg(P)
    ok = true; msg = "";

    % Nonnegative fields (existing four)
    nonnegFields = {'threshold_brain','threshold_scalp','sigma_brain','sigma_scalp'};
    for k=1:numel(nonnegFields)
        fld = nonnegFields{k};
        v = P.(fld);
        if ~isscalar(v) || ~isfinite(v) || v < 0
            ok=false; msg=sprintf('%s must be a finite number >= 0.', fld); return;
        end
    end

    % threshold_sensitivity: any finite scalar allowed (can be negative)
    v = P.threshold_sensitivity;
    if ~isscalar(v) || ~isfinite(v)
        ok=false; msg='threshold_sensitivity must be a finite number.'; return;
    end
end
