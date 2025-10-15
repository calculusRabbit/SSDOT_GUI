function plotRaw(uiAx, runObj, userSel)
    acq = runObj.acquired;

    % Get time vector and data matrix
    y = acq.GetDataTimeSeries();
    t = acq.GetTime();


    % Measurement list (numeric matrix)
    ml = acq.GetMeasList();   % columns: [src det type wlIdx dataTypeIdx ...]

    % Extract user selection
    src = userSel.source;
    det = userSel.detector;


    for i = 1:numel(userSel.selectSignal)
        wlIdx = userSel.selectSignal(i);

        % Find matching measurement column
        col = find(ml(:,1) == src & ml(:,2) == det & ml(:,4) == wlIdx, 1, 'first');

        % Plot the raw
        label = sprintf('Run %d - Tx%d-Rx%d - Raw - %s', ...
            runObj.iRun, src, det, userSel.waveLengthLabel{wlIdx});
        plot(uiAx, t, y(:,col), 'LineWidth', 1.3, 'DisplayName', label);
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
