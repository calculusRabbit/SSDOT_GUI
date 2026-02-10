function plotOD(uiAx, runObj, userSel)
    source = userSel.source;
    detector = userSel.detector;
    % get optical density by using homer3 GetVar
    dod = runObj.procStream.output.GetVar('dod'); % GetVar using homer3 function
    % channel info (Tx,Rx,wavelength)
    ml = dod.measurementList;

    % build arrays of Tx, Rx, and wavelength index
    % To search which column in the DataTimeSeries matix matches the user’s chosen (Tx, Rx, wavelength).
    sourceArrayOD = zeros(numel(ml), 1);
    detectorArrayOD = zeros(numel(ml), 1);
    waveLengthArrayOD = zeros(numel(ml), 1);

    
    %fill the arrays
    for i = 1:(numel(ml))
        sourceArrayOD(i) = ml(i).sourceIndex;
        detectorArrayOD(i) = ml(i).detectorIndex;
        waveLengthArrayOD(i) = ml(i).wavelengthIndex;
    end

    x = dod.GetTime(); % time vector
    y = dod.GetDataTimeSeries(); % OD data matrix

    % loop over user selected wavelength
    for i = 1:numel(userSel.selectSignal)
        % eg: map 850nm -> index        
        % find the column index that match
        % source, detector, wavelength
        col = find(source == sourceArrayOD & ...
            detector == detectorArrayOD &  ...
            userSel.selectSignal(i) == waveLengthArrayOD, 1, 'first');

        % plot OD curve
        dataTypeLabel = ml(col).dataTypeLabel;
        name = sprintf('Run %d — Sx%d–Dx%d — %s - %s', runObj.iRun, source, detector, dataTypeLabel, userSel.waveLengthLabel{userSel.selectSignal(i)});
        plot(uiAx, x, y(:,col), 'LineWidth', 1.3, 'DisplayName', name);
    end

    %% debug:
    if (runObj.iRun == 2)
        return;
    end


    % Overlay stim vertical lines
    stimList = runObj.GetStim();

    %colorOrder = get(uiAx, 'ColorOrder');

    for i = 1:numel(stimList)
        data = stimList(i).GetData();
        if isempty(data)
            continue
        end

        % first column is the onset timeline
        onsetTime = data(:,1);
        
        % Generate unique color based on run and condition
        hue = mod((runObj.iRun - 1) * numel(stimList) + i - 1, 12) / 12;
        color = hsv2rgb([hue, 0.8, 0.9]); 

        for j = 1:numel(onsetTime)
            if j == 1
            xline(uiAx, onsetTime(j), '--', ...
                'LineWidth', 1.3, ...
                'Color', color, ...
                'DisplayName', sprintf('Run %d: %s', runObj.iRun, stimList(i).name));
            else
            xline(uiAx, onsetTime(j), '--', ...
                'LineWidth', 1.3, ...
                'Color', color, ...
                'HandleVisibility', 'off');
            end
        end
    end
end
