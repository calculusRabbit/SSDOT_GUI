function plotODAvg(uiAx, obj, userSel)
    dodAvg = obj.procStream.output.GetVar('dodAvg');
    if isempty(dodAvg)
        warning('plotODAvg: dodAvg not found in output.');
        fig = ancestor(uiAx, 'figure');   
        uialert(fig, 'plotODAvg: dodAvg not found in output.', 'Warning');
        return;
    end

    x = dod.GetTime(); % time vector
    y = dod.GetDataTimeSeries(); % OD data matrix

    ml = dodAvg.measurementList;
    for i = 1:numel(userSel.selectSignal)
        for j = 1:numel(ml)
            condName = userSel.condName;
            condLabel = obj.CondNames{ml(j).dataTypeIndex};
            if (ml(j).sourceIndex == userSel.source &&...
                ml(j).detectorIndex == userSel.detector &&...
                ml(j).wavelengthIndex == userSel.selectSignal(j) && ...
                strcmp(condName, condLabel))
                col = j;
                break;
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

            dataTypeLabel = ml(col).dataTypeLabel;
            name = sprintf('%s %d — Sx%d-Dx%d — %s - %s',obj.type, idx, userSel.source, userSel.detector, dataTypeLabel, userSel.waveLengthLabel{userSel.selectSignal(i)});
            plot(uiAx, x, y(:,col), 'LineWidth', 1.3, 'DisplayName', name);
        end
    end
end