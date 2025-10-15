function plotConcAvg(uiAx, obj, userSel)
    dcAvg = obj.procStream.output.GetVar('dcAvg');
    if isempty(dcAvg)
        warning('plotdcAvg: dcAvg not found in output.');
        fig = ancestor(uiAx, 'figure');   
        uialert(fig, 'plotdcAvg: dcAvg not found in output.', 'Warning');
        return;
    end

    ml = dcAvg.measurementList;

    for i = 1:numel(ml)
        condName = userSel.condName;
        condLabel = obj.CondNames{ml(i).dataTypeIndex};
        if (userSel.source == ml(i).sourceIndex && userSel.detector == ml(i).detectorIndex && (strcmp(condName, condLabel)));
            hbo_idx = i;
            break;
        end
    end

    if strcmp('run', obj.type)
        idx = obj.iRun;
    elseif strcmp('sess', obj.type)
        idx = obj.iSess;
    elseif strcmp('subj', obj.type)
        idx = obj.iSubj;
    else
        idx = obj.iGroup;
    end

    x = dcAvg.GetTime();
    y = dcAvg.GetDataTimeSeries();
    

    for i = 1:numel(userSel.selectSignal)
        k = userSel.selectSignal(i);

        dataTypeLabel = dcAvg.measurementList(hbo_idx + k).dataTypeLabel;
        name = sprintf('%s %d — %s — Tx%d-Rx%d', obj.type, idx, dataTypeLabel, userSel.source, userSel.detector);
        plot(uiAx, x, y(:, hbo_idx + k), 'LineWidth', 1.3, 'DisplayName', name);
    end
end