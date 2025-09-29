function plotConc(uiAx, runObj, userSel)
    runObj.Load();
    dc = runObj.procStream.output.GetVar('dc');
    ml = dc.measurementList;

    % find HbO index
    for i = 1:numel(ml)
        if (userSel.source == ml(i).sourceIndex && userSel.detector == ml(i).detectorIndex)
            hbo_idx = i; 
            break;
        end
    end

    x = dc.GetTime();
    y = dc.GetDataTimeSeries();

    for i = 1:numel(userSel.selectSignal)
        k = userSel.selectSignal(i);

        dataTypeLabel = dc.measurementList(hbo_idx + k).dataTypeLabel;
        name = sprintf('Run %d — %s — Tx%d–Rx%d', runObj.iRun, dataTypeLabel, userSel.source, userSel.detector);
        plot(uiAx, x, y(:, hbo_idx + k), 'LineWidth', 1.3, 'DisplayName', name);
    end
end
