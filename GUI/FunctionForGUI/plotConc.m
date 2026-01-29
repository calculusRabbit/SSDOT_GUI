function plotConc(uiAx, runObj, userSel)
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
        name = sprintf('Run %d — %s — Sx%d–Dx%d', runObj.iRun, dataTypeLabel, userSel.source, userSel.detector);
        plot(uiAx, x, y(:, hbo_idx + k), 'LineWidth', 1.3, 'DisplayName', name);
    end


    % Overlay stim vertical lines
    stimList = runObj.GetStim();

    colorOrder = get(uiAx, 'ColorOrder');

    for i = 1:numel(stimList)
        data = stimList(i).GetData();
        if isempty(data)
            continue
        end

        % first column is the onset timeline
        onsetTime = data(:,1);
        color = colorOrder(i, :);

        for j = 1:numel(onsetTime)
            if j == 1
            xline(uiAx, onsetTime(j), '--', ...
                'LineWidth', 1.3, ...
                'Color', color, ...
                'DisplayName', stimList(i).name);
            else
            xline(uiAx, onsetTime(j), '--', ...
                'LineWidth', 1.3, ...
                'Color', color, ...
                'HandleVisibility', 'off');
            end
        end
    end
end
