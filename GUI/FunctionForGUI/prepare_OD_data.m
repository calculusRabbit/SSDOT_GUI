function prepare_OD_data(app, ssdot)
    %% Calculate Y
    % Sensitivity_Matrix = get_global_variable(app, 'Sensitivity_Matrix');
    Sensitivity_Matrix = ssdot.getVar('Sensitivity_Matrix');
    channels = Sensitivity_Matrix.channels;
    % VU this needs to look at which run is selected
    % data = load('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\sess01\2025-08-29_001.mat');
    % OD  = data.output.dod.dataTimeSeries;
    data = ssdot.getVar('RunSelected').procStream;
    OD = data.output.dod.dataTimeSeries;

    [time_points,n_channels] = size(OD);
    OD_w1 = OD(:,channels)';
    OD_w2 = OD(:,channels + size(OD,2)/2)';
    Y_w1 = reshape(OD_w1,[],1);
    Y_w2 = reshape(OD_w2,[],1);
    Y = [Y_w1;Y_w2];
    % Y_w1 t1 ch1
    %         ch2
    %          :
    %         chn
    %      t2 ch1
    %         ch2
    %          :
    %         chn
    % Y_w2
    % calculated_mat.Y = Y;
    % save_global_data(app,'Y',Y); % Vu you need to change here
    ssdot.setVar('Y', Y);

    fprintf('finish Y\n');

    active_channel_number = length(channels)*2;
    %% Calculate SS OD time series
    % VU: where to find
    % C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\ProcStreamFunctionsSummary.txt
    % [dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, [-2.0 18.0], 1, 2, [0.1 3.0 7.0 0.5 3.0 7.0], 22.0, 2 (should be the 2 here), 3, 0);
    flagNuisanceRMethod = 2;
    rhoSD_ssThresh = 22;

    if flagNuisanceRMethod == 2
        % atlasViewer = get_global_variable(app, 'AtlasState');
        atlasViewer = ssdot.getVar('AtlasState');
        probe = atlasViewer.probe;
        SD = convertProbe2SD(probe);

        ml = SD.MeasList;
        lst = find(ml(:,4)==1 & SD.MeasListAct==1);
        rhoSD = zeros(length(lst),1);
        posM = zeros(length(lst),3);
        for iML = 1:length(lst)
            rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
            posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
        end
        shortSepChLst = lst(find(rhoSD<rhoSD_ssThresh));
        % Vu: how to get mlActAuto:
        % C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\sess01\2025-08-29_001.mat
        output = load('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\sess01\2025-08-29_001.mat','output');
        mlActAuto = output.output.misc.mlActAuto;
        if size(ml,1) ~= size( mlActAuto{1,1},1)
            mlActAuto{1,1}([51 102],:) = [];
        end
        activeChLst = find(ml(:,4)==1 & mlActAuto{1,1}==1);
        [channels_short,~]=intersect(activeChLst,shortSepChLst); % channels_short is active short channels

        OD_SS_w1 = OD(:, channels_short);
        OD_SS_w2 = OD(:, n_channels/2 + channels_short);
        OD_SS_w1_avg = mean(OD_SS_w1,2);
        OD_SS_w2_avg = mean(OD_SS_w2,2);

        
        OD_SS = zeros(time_points*active_channel_number,active_channel_number);
        for i = 1:active_channel_number/2
            OD_SS(i:active_channel_number/2:time_points*active_channel_number/2,i) = OD_SS_w1_avg;
        end
        for i = 1:active_channel_number/2
            OD_SS(time_points*active_channel_number/2 + i:active_channel_number/2:time_points*active_channel_number, active_channel_number/2 + i) = OD_SS_w2_avg;
        end
    %% currently only maintain when we average the shorter channels as our SS, in the future can expand
    end
    %% Calculate drift OD_drift
    % VU: where to find
    % C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\ProcStreamFunctionsSummary.txt
    % [dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, [-2.0 18.0], 1, 2, [0.1 3.0 7.0 0.5 3.0 7.0], 22.0, 2 , 3 (should be the 3 here), 0);
    driftOrder = 3;
    % T = get_global_variable(app,'T');
    T = ssdot.getVar('T');
    D = make_D(T, driftOrder);
    [time_points,n_drifter] = size(D);
    OD_drift = zeros(active_channel_number*time_points,active_channel_number*n_drifter);
    for i = 1:active_channel_number/2
        for j = 1:n_drifter
            OD_drift(i:active_channel_number/2:time_points*active_channel_number/2,(i-1)*n_drifter+j) = D(:,j);
        end
    end
    for i = 1:active_channel_number/2
        for j = 1:n_drifter
            OD_drift(time_points*active_channel_number/2 + i:active_channel_number/2:time_points*active_channel_number, active_channel_number*n_drifter/2 + (i-1)*n_drifter+j) = D(:,j);
        end
    end
    % save_global_data(app, 'OD_SS', OD_SS)
    ssdot.setVar('OD_SS', OD_SS);
    % save_global_data(app, 'OD_drift', OD_drift)
    ssdot.setVar('OD_drift');
end