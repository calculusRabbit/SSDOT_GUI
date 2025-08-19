function calculated_mat = Proc_stream_v7(snirf_file, mat_file, cfg, calculated_mat)

if isempty(calculated_mat)
    mask_threshold = cfg.mask_threshold;
    threshold_brain = cfg.threshold_brain;
    threshold_scalp = cfg.threshold_scalp;
    sigma_brain = cfg.sigma_brain;
    sigma_scalp = cfg.sigma_scalp;
    rhoSD_ssThresh = cfg.rhoSD_ssThresh; % 7 mm
    device = cfg.device; % device could be 'gpu' or 'cpu'
    
    Proc_data = process_in_homer_fromOD(snirf_file, mat_file, cfg);
    calculated_mat.Proc_data = Proc_data;
    
    %% Calculate A
    
    Sensitivity_Matrix= Get_A_dot(cfg.av_data_path, rhoSD_ssThresh.av, Proc_data.mlActAuto, cfg.spatially_regu);
    calculated_mat.Adot = Sensitivity_Matrix.Adot;
    calculated_mat.Adot_scalp = Sensitivity_Matrix.Adot_scalp;
    calculated_mat.E = Sensitivity_Matrix.E;
    calculated_mat.channels = Sensitivity_Matrix.channels;
    calculated_mat.shortSepChLst = Sensitivity_Matrix.shortSepChLst;
    if cfg.spatially_regu
        calculated_mat.ll = Sensitivity_Matrix.ll;
        calculated_mat.ll_scalp = Sensitivity_Matrix.ll_scalp;
    end
    At_file = 'atlasViewer.mat';

    M = Make_mask(mask_threshold, Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);
    calculated_mat.M = M;
    % here we use the spatially regularized A matrix 
    A = Make_A_matrix(Sensitivity_Matrix.Adot, Sensitivity_Matrix.Adot_scalp, Sensitivity_Matrix.E,  M);
    calculated_mat.A = A;

    %% Calculate T
    T = Make_T_matrix(Proc_data);
    calculated_mat.T = T;
    
    %% Calculate G
    G = Make_G_matrix([cfg.av_data_path,At_file], M, threshold_brain, threshold_scalp, sigma_brain, sigma_scalp, Sensitivity_Matrix.Adot_orig, Sensitivity_Matrix.Adot_scalp_orig);%% David Beeler
    calculated_mat.G = G;

    %% Calculate Y 
    channels = Sensitivity_Matrix.channels;
    OD = Proc_data.dod.dataTimeSeries;
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
    calculated_mat.Y = Y;
    active_channel_number = length(channels)*2;
    %% Calculate SS OD time series
    if cfg.SS_flag == 1
        At_file = 'atlasViewer.mat';
        atlasViewer = load([cfg.av_data_path,At_file]);
        warning off
        
        probe = atlasViewer.probe;
        SD = convertProbe2SD(probe);

        ml          = SD.MeasList;
        lst = find(ml(:,4)==1 & SD.MeasListAct==1);
        rhoSD = zeros(length(lst),1);
        posM = zeros(length(lst),3);
        for iML = 1:length(lst)
            rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
            posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
        end
        shortSepChLst = lst(find(rhoSD<rhoSD_ssThresh.homer));
        ml = SD.MeasList;
        mlActAuto = Proc_data.mlActAuto;
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
    elseif cfg.SS_flag == 0
        OD_SS = [];
    elseif cfg.SS_flag == 2
        load(mat_file, 'OD');
        OD_SS_w1_avg = OD(:,51);
        OD_SS_w2_avg = OD(:,102);
        
        OD_SS = zeros(time_points*active_channel_number,active_channel_number);
        for i = 1:active_channel_number/2
            OD_SS(i:active_channel_number/2:time_points*active_channel_number/2,i) = OD_SS_w1_avg;
        end
        for i = 1:active_channel_number/2
            OD_SS(time_points*active_channel_number/2 + i:active_channel_number/2:time_points*active_channel_number, active_channel_number/2 + i) = OD_SS_w2_avg;
        end
    end
    %% Calculate drift OD_drift
    if cfg.drift_flag == 1
        if isfield(cfg,'driftOrder')
            driftOrder = cfg.driftOrder;
        else
            driftOrder = 0;
        end
        D = make_D(calculated_mat.T, driftOrder);
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
    elseif cfg.drift_flag == 0
        n_drifter = 0;
        OD_drift = [];
    end
    
    %% Calculate H: H = [OD_SS, OD_drift, H_brain, H_scalp]
    [H_brain,H_scalp] = make_H_matrix(A, G, T, Y, device, cfg);
    calculated_mat.H_scalp = H_scalp;
    calculated_mat.H_brain = H_brain;
    calculated_mat.OD_SS = OD_SS;
    calculated_mat.OD_drift = OD_drift;
    calculated_mat.n_channels = n_channels;
    calculated_mat.n_drifter = n_drifter;
    %% calculate HTH and HTY
    chunk = 1024;
    if cfg.regularization == 2
        Observations = [OD_w1', OD_w2'];
        R = estimateMeasureNoise(Observations);
        R = repmat(R, 1, size(OD_w1,2));
        R = spdiags(R',0,speye(length(R),length(R)));
        R = inv(R);
    else
        R = [];
    end
    if cfg.regularization ~= 3
        [HTY, HTH] = make_HTY_and_HTH(H_brain,    H_scalp,    OD_SS, OD_drift, Y, chunk, device, R);
        calculated_mat.HTH = HTH;
        calculated_mat.HTY = HTY;
    end
end
%% inverse

if cfg.regularization == 3
    H = [H_brain, H_scalp, OD_SS, OD_drift]; 
    [u,s,v]=svds(H,size(H,2)); 
    max_sing = max(s(:));
    alpha = cfg.alpha * max_sing;
    H_new = u * sqrtm(s*s + alpha^2*eye(size(s))) *v';
    b = inv(H_new'*H_new)* (H_new' * Y) ;
    calculated_mat.b = b;
else
    HTH = calculated_mat.HTH;
    HTY = calculated_mat.HTY;
    if ~isfield(cfg,'beta')
        % first regularization method
        alpha = cfg.alpha;
        n_HTH = size(HTH,2);

        r_alpha = diag( alpha * eigs(HTH,1).* ones(1,n_HTH));
        b = HTY'/(HTH + r_alpha);

        calculated_mat.b = b;
    else
        % second regularization method
        alpha = cfg.alpha;
        beta = cfg.beta;
        H_brain = calculated_mat.H_brain;
        H_scalp = calculated_mat.H_scalp;

        n_brain = size(H_brain,2);
        n_scalp = size(H_scalp,2);
        HTH_brain = HTH(1:n_brain + n_scalp, 1:n_brain + n_scalp);
        HTH_scalp = HTH(n_brain + n_scalp + 1:end,n_brain + n_scalp + 1:end);
        n_brain = size(HTH_brain,2);
        n_others = size(HTH_scalp,2);
        r_alpha = diag([alpha * eigs(HTH_brain,1).* ones(1,n_brain), zeros(1,n_others)]);
        if isempty(HTH_scalp)
            r_beta = 0;
        else
            r_beta = diag([zeros(1,n_brain), beta * eigs(HTH_scalp,1).* ones(1,n_others)]);
        end
        b = HTY'/(HTH + r_beta + r_alpha);
        calculated_mat.b = b;
    end
end

