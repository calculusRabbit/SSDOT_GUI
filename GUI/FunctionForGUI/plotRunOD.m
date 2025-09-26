function visualizeRunConc(uiAx, runObj, source, detector, selectSignal, runIdx)
    runObj.Load();
    dc = runObj.procStream.output.GetVar('dc');
    ml = dc.measurementList;

    % find HbO index
    for i = 1:numel(ml)
        if source == ml(i).sourceIndex && detector == ml(i).detectorIndex
            hbo_idx = i; break;
        end
    end

    x = dc.GetTime();
    y = dc.GetDataTimeSeries();

    for kSel = 1:numel(selectSignal)
        sel = lower(selectSignal{kSel});
        if strcmp(sel,'hbo'), k = 0;
        elseif strcmp(sel,'hbr'), k = 1;
        elseif strcmp(sel,'hbt'), k = 2;
        else, continue;
        end

        tag = dc.measurementList(hbo_idx + k).dataTypeLabel;
        name = sprintf('Run %d — %s — Tx%d–Rx%d', runIdx, tag, source, detector);
        plot(uiAx, x, y(:, hbo_idx + k), 'LineWidth', 1.3, 'DisplayName', name);
    end
end
