function plot_single_time_base()
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

    function get_x1()
        % where you can find x1 from selected run
        % check if a run is selected - list
        % check if run level is selected - ratio buttone
        run_number = get_selected_run();

    end

    % check if we are on the run level;

    idxBasis = get_idxBasis();
    paramsBasis = get_paramsBasis();
    [tHRF,dt] = get_tHRF();
    ntHRF = length(tHRF);
    trange = get_trange();
    tbasis = construct_basis_function(idxBasis,paramsBasis,trange,ntHRF, tHRF,dt);
    % need to do convolution
    x1 = get_x1();
    conv(x1 , x2)
end