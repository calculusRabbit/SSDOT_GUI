function plotOD(uiAx, runObj, userSel)

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
        waveLengthIdx = app.nmLabelToIndex(userSel.selectSignal{i});
        
        % find the column index that match
        % source, detector, wavelength
        col = find(source == sourceArrayOD & ...
            detector == detectorArrayOD &  ...
            waveLengthIdx == waveLengthArrayOD, 1, 'first');

        % plot OD curve
        name = sprintf('Run %d — %s — Tx%d–Rx%d', runObj.iRun, dataTypeLabel, source, detector);
        plot(uiAx, x, y(:,col), 'LineWidth', 1.3, 'DisplayName', name);
    end
end
