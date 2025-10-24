function [cfgOut, accepted] = userOptionGUI(cfgIn)
% Modal dialog to edit processing thresholds/sigmas.
% Returns updated cfg and accepted=true if OK pressed.

    if nargin<1 || isempty(cfgIn)
        cfgIn = struct('threshold_brain',5,'threshold_scalp',5,'sigma_brain',5,'sigma_scalp',5);
    end
    cfg = cfgIn; 
    accepted = false;

    fig = uifigure('Name','UserOptionGUI', ...
                   'Position',[100 100 460 210], ...
                   'WindowStyle','modal', ...
                   'CloseRequestFcn', @(source,event) doClose());

    g = uigridlayout(fig,[5 2]);
    g.RowHeight   = {30,30,30,30,44};
    g.ColumnWidth = {200,'1x'};
    g.Padding     = [12 12 12 12];
    g.RowSpacing  = 8; g.ColumnSpacing = 10;

    % ---- labels + fields ----
    uilabel(g,'Text','threshold_brain','HorizontalAlignment','right');
    efTB = uieditfield(g,'numeric','Value',cfg.threshold_brain,'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','threshold_scalp','HorizontalAlignment','right');
    efTS = uieditfield(g,'numeric','Value',cfg.threshold_scalp,'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','sigma_brain','HorizontalAlignment','right');
    efSB = uieditfield(g,'numeric','Value',cfg.sigma_brain,'LowerLimitInclusive','on','LowerLimit',0);

    uilabel(g,'Text','sigma_scalp','HorizontalAlignment','right');
    efSS = uieditfield(g,'numeric','Value',cfg.sigma_scalp,'LowerLimitInclusive','on','LowerLimit',0);

    % ---- buttons ----
    uibutton(g,'Text','Cancel','ButtonPushedFcn',@(s,e) doClose(), ...
        'Layout',struct('Row',5,'Column',1));
    uibutton(g,'Text','OK','ButtonPushedFcn',@(s,e) onOK(), ...
        'Layout',struct('Row',5,'Column',2));

    uiwait(fig);

    if accepted, cfgOut = cfg; else, cfgOut = cfgIn; end

    % ===== nested helpers =====
    function onOK()
        % read
        P.threshold_brain = efTB.Value;
        P.threshold_scalp = efTS.Value;
        P.sigma_brain     = efSB.Value;
        P.sigma_scalp     = efSS.Value;

        % validate (adjust rules as you like)
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
    fields = {'threshold_brain','threshold_scalp','sigma_brain','sigma_scalp'};
    for k=1:numel(fields)
        v = P.(fields{k});
        if ~isscalar(v) || ~isfinite(v) || v < 0
            ok=false; msg=sprintf('%s must be a finite number >= 0.', fields{k}); return;
        end
    end

end
