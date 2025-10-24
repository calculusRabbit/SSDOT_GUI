function plot_single_time_base(app, ssdot)
    function idxBasis = get_idxBasis()
        % Vu: in the derivatives/homer/ProcStreamFunctionsSummary.txt you can
        % find idx Basis in line 36: 
        % [dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = 
        % hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, 
        % [-2.0 18.0], 1, 2, [0.1 3.0 7.0 0.5 3.0 7.0], 22.0, 2, 3, 0);

        idxBasis = 2; %third one
    end
    
    function paramsBasis = get_paramsBasis()
        % Vu: in the derivatives/homer/ProcStreamFunctionsSummary.txt you can
        % find idx Basis in line 36: 
        % [dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = 
        % hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, 
        % [-2.0 18.0], 1, 2, [0.1 3.0 7.0 0.5 3.0 7.0], 22.0, 2, 3, 0);
        paramsBasis = [0.1 3.0 7.0 0.5 3.0 7.0];% forth one
    end
    
    function [tHRF,dt] = get_tHRF()
        % Vu: please revise here so that it will loop the files in sess01
        data =  load('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\sess01\2025-08-29_001.mat','output');
        
        
        output = data.output;
        t = output.dod.time;

        dt = t(2) - t(1);
        trange = get_trange();
        nPre = round(trange(1)/dt);
        nPost = round(trange(2)/dt);
        tHRF = (1*nPre*dt:dt:nPost*dt)';
    end
    
    function trange = get_trange()
        % Vu: in the derivatives/homer/ProcStreamFunctionsSummary.txt you can
        % find idx Basis in line 36: 
        % [dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = 
        % hmrR_GLM(dc, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, 
        % [-2.0 18.0], 1, 2, [0.1 3.0 7.0 0.5 3.0 7.0], 22.0, 2, 3, 0);
        trange = [-2.0 18.0]; % first one
    end

    function nT = get_nT()
        % where you can find nT from selected run: derivatives folder
        % check if a run is selected - list
        % check if run level is selected - ratio buttone
        data =  load('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\derivatives\homer\sess01\2025-08-29_001.mat','output');
        output = data.output;
        t = output.dod.time;
        nT = length(t);
    end

    function nCond = get_nCond()
        % where you can find snirf from selected run
        snirf = SnirfClass('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\sess01\2025-08-29_001.snirf');
        t = snirf.GetTimeCombined();
        stimStates = snirf.GetStims(t);
        lstCond = find(sum(stimStates == 1, 1) > 0);
        nCond = length(lstCond); % size(stimStates,2);
    end
    function onset = get_onset()
        
        %%%%%%%%%%%%%%%%
        % Prune good stim, generate onset matrix
        %%%%%%%%%%%%%%%%
        % Get only indices of conditions with any stimStates that are 1
        % where you can find snirf from selected run
        snirf = SnirfClass('C:\Users\q582z568\Documents\SSDOT_GUI\ST001\sess01\2025-08-29_001.snirf');
        t = snirf.GetTimeCombined();
        stimStates = snirf.GetStims(t);
        dt = t(2) - t(1);
        trange = get_trange();
        nPre = round(trange(1)/dt);
        nPost = round(trange(2)/dt);
        nTpts = length(snirf.data.time);
        stim = snirf.stim;
        
        iBlk = 1;
        lstCond = find(sum(stimStates == 1, 1) > 0);
        nCond = length(lstCond); % size(stimStates,2);
        nTrials{iBlk} = zeros(nCond,1);
        onset = zeros(nT, nCond);
        avg_pulses = {};
        
        for iCond = 1:nCond
            lstT = find(stimStates(:, lstCond(iCond)) == 1);  % Indices of stims enabled (== 1)
            lstp = find((lstT+nPre) >= 1 & (lstT+nPost) <= nTpts);  % Indices of stims not clipped by signal
            lst = lstT(lstp);
            nTrials{iBlk}(iCond) = length(lst);
            % Generate basis boxcars of stim amplitude and duration
            starts = lst+nPre;
            if ~isempty(stim(lstCond(iCond)).data)
                durations = stim(lstCond(iCond)).data(:, 2);
                amplitudes = stim(lstCond(iCond)).data(:, 3);
                avg_pulses{iCond} = ones(round(mean(durations) / dt), 1); %#ok<AGROW>
                for i = 1:length(starts)
                    if idxBasis == 1  % Gaussian has no duration T (yet)
                       pulse_duration = 1; 
                    else
                       pulse_duration = round(durations(i) / dt); 
                    end
                    pulse = (amplitudes(i) / pulse_duration) * ones(pulse_duration, 1);
                    onset(starts(i):starts(i) + pulse_duration - 1, iCond) = onset(starts(i):starts(i) + pulse_duration - 1, iCond) + pulse;
                end
            end
        end
    end



    idxBasis = get_idxBasis();
    paramsBasis = get_paramsBasis();
    [tHRF,dt] = get_tHRF();
    ntHRF = length(tHRF);
    trange = get_trange();
    [tbasis, nB] = construct_basis_function(idxBasis,paramsBasis,trange,ntHRF, tHRF,dt);
     
    figure; plot(squeeze(tbasis(:,:,1)));title("single HbO tbasis")
    figure; plot(squeeze(tbasis(:,:,2)));title("HbR tbasis")
    
    nT = get_nT();
    nCond = get_nCond();
    onset = get_onset();

    dA = construct_design_matrix(nT,nB,nCond,onset,tbasis);
    T = Make_T_matrix(dA,tbasis);
    fprintf('finish T\n')
    figure; plot(squeeze(T.T_HbO_brain));title("convolved HbO tbasis"); legend('cond1','cond2');
    figure; plot(squeeze(T.T_HbO_brain));title("convolved HbR tbasis"); legend('cond1','cond2');
    % VU help me save T to global space
    % save_global_data(app, 'T', T)
    ssdot.setVar('T', T);
end