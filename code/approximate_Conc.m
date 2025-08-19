function Conc = approximate_Conc(calculated_mat, cfg)

canUseGPU=true;

if strcmp(cfg.device,'gpu')
    try
%         fprintf('Finding GPU...\n')
        if gpuDeviceCount > 1
            for ii = 1:gpuDeviceCount
                g = gpuDevice(ii);
                fprintf(1,'Device %i has ComputeCapability %s \n', ...
                        g.Index,g.ComputeCapability)
            end
            prompt = 'Which GPU? ';
            x = input(prompt);
            d = gpuDevice(x);
        else
            d = gpuDevice;
        end
    catch
        canUseGPU=false;
        fprintf('no GPU. Using cpu now\n')
    end
    n_batch = cfg.n_batch_brain;
    H_finished = false;
    i = 1;
    while i <= 100
        fprintf('app Conc batch %d\n', n_batch)
        try
            Conc = calculate_Conc(calculated_mat.Proc_data, calculated_mat.G, calculated_mat.T, calculated_mat.b, n_batch, cfg);
            H_finished = true;
            break
        catch
            n_batch = n_batch*2;
            i = i + 1;
            continue
        end

    end
    if H_finished == false
        fprintf('Cannot fit into GPU memory even with %d batches\n', n_batch)
        return
    end
    
elseif strcmp(cfg.device, 'cpu') || ~canUseGPU
    cond_num = length(calculated_mat.Proc_data.stim);
    for i_cond = 1:cond_num
        if cfg.with_scalp == 1  
            T_HbO_brain = calculated_mat.T.T_HbO_brain;
            T_HbR_brain = calculated_mat.T.T_HbR_brain;
        else
            T_HbO_brain = calculated_mat.T.t_HbO_brain;
            T_HbR_brain = calculated_mat.T.t_HbR_brain;
        end
    
        n_g = size(calculated_mat.G.G_brain,1);
        tbasis_num = size(T_HbO_brain,2);
        n_t = tbasis_num * cond_num;
        index = 1:n_g*n_t;
        index = reshape(index, cond_num*tbasis_num, n_g);
        index = index((i_cond-1)*tbasis_num+1:i_cond*tbasis_num,:);
        index = reshape(index, [],1);
        beta_brain_HbO = calculated_mat.b(index); % H_brain [g1[c1[t1-tn] c2],g2[c1[t1-tn] c2],...]
        beta_brain_HbR = calculated_mat.b(n_g*n_t + index);

        G_brain = calculated_mat.G.G_brain;

        Conc_HbO = zeros(size(calculated_mat.G.G_brain,2),size(T_HbO_brain,1));
        Conc_HbR = zeros(size(calculated_mat.G.G_brain,2),size(T_HbR_brain,1));

        j = 1;

        for i_g = 1:n_g
            for i_t = 1:n_t
                gt_HbO = kron(G_brain(i_g,:),T_HbO_brain(:,i_t));
                gt_HbR = kron(G_brain(i_g,:),T_HbR_brain(:,i_t));

                Conc_HbO = Conc_HbO + gt_HbO'.*beta_brain_HbO(j);
                Conc_HbR = Conc_HbR + gt_HbR'.*beta_brain_HbR(j);
                j = j + 1;
            end
        end

        Conc.intensity_HbO{i_cond} = Conc_HbO;
        Conc.intensity_HbR{i_cond} = Conc_HbR;
    
        if cfg.with_scalp == 1
            n_g_scalp = size(calculated_mat.G.G_scalp,1);
            tbasis_num = size(calculated_mat.T.T_HbO_scalp,2);
            n_t_scalp  = tbasis_num*cond_num;

            index = 1:n_g_scalp*n_t_scalp;
            index = reshape(index, cond_num*tbasis_num, n_g_scalp);
            index = index((i_cond-1)*tbasis_num+1:i_cond*tbasis_num,:);
            index = reshape(index, [],1);
            beta_scalp_HbO = calculated_mat.b(2*n_g*n_t + index); % H_scalp [g1[c1[t1-tn] c2],g2[c1[t1-tn] c2],...]
            beta_scalp_HbR = calculated_mat.b(2*n_g*n_t + n_g_scalp*n_t_scalp + index);

            Conc_HbO_scalp = zeros(size(calculated_mat.G.G_scalp,2),size(calculated_mat.T.T_HbO_scalp,1));
            Conc_HbR_scalp = zeros(size(calculated_mat.G.G_scalp,2),size(calculated_mat.T.T_HbR_scalp,1));

            G_scalp = calculated_mat.G.G_scalp;
            T_HbO_scalp = calculated_mat.T.T_HbO_scalp;
            T_HbR_scalp = calculated_mat.T.T_HbR_scalp;

            j = 1;

            for i_g = 1:n_g_scalp
                for i_t = 1:n_t_scalp
                    gt_HbO = kron(G_scalp(i_g,:),T_HbO_scalp(:,i_t));
                    gt_HbR = kron(G_scalp(i_g,:),T_HbR_scalp(:,i_t));

                    Conc_HbO_scalp = Conc_HbO_scalp + gt_HbO'.*beta_scalp_HbO(j);
                    Conc_HbR_scalp = Conc_HbR_scalp + gt_HbR'.*beta_scalp_HbR(j);
                    j = j + 1;
                end
            end

            Conc.intensity_HbO_scalp{i_cond} = Conc_HbO_scalp;
            Conc.intensity_HbR_scalp{i_cond} = Conc_HbR_scalp;
        end
    end
end
end

function Conc = calculate_Conc(Proc_data,G, T, beta, n_batch,cfg)
cond_num = length(Proc_data.stim);

% [time_points, [Cond1[tB1...tBx], Cond2[tB1...tBx]...], Concentration]
for i_cond = 1:cond_num
    if cfg.with_scalp == 1 || (isfield(cfg,'whole_time') && cfg.whole_time == 1)
        T_HbO_brain = T.T_HbO_brain;
        T_HbR_brain = T.T_HbR_brain;
    else
        T_HbO_brain = T.t_HbO_brain;
        T_HbR_brain = T.t_HbR_brain;
    end
    n_g = size(G.G_brain,1);
    tbasis_num = size(T.t_HbO_brain,2);
    n_t = tbasis_num*cond_num;
    index = 1:n_g*n_t;
    index = reshape(index, cond_num*tbasis_num, n_g);
    index = index((i_cond-1)*tbasis_num+1:i_cond*tbasis_num,:);
    index = reshape(index, [],1);
    beta_brain_HbO = beta(index); % H_brain [g1[c1[t1-tn] c2],g2[c1[t1-tn] c2],...]
    beta_brain_HbR = beta(n_g*n_t + index);
    
    Conc_HbO = zeros(size(G.G_brain,2),size(T_HbO_brain,1),'gpuArray');
    Conc_HbR = zeros(size(G.G_brain,2),size(T_HbO_brain,1),'gpuArray');

    n_g_batch = floor(n_g/n_batch);

    i_g = 1;
    while i_g <= n_g
        g_start     =   i_g;
        g_end       =   i_g + n_g_batch - 1 ;
        if g_end > n_g
            g_end = n_g;
        end
        g_index     =   g_start:g_end;

        Gm = gpuArray(G.G_brain(g_index,:));
        T_HbO = gpuArray(T_HbO_brain);
        T_HbR = gpuArray(T_HbR_brain);
        beta_HbO = beta_brain_HbO((g_start-1)*tbasis_num+1 : g_end*tbasis_num);
        beta_HbR = beta_brain_HbR((g_start-1)*tbasis_num+1 : g_end*tbasis_num);

        HbO = kron(reshape(Gm',[],1),reshape(T_HbO,1,[]));
        HbR = kron(reshape(Gm',[],1),reshape(T_HbR,1,[]));
        [m, n] = size(Gm);
        [k, l] = size(T_HbO);
        h = 1;
        for j = 1:m
            for i = 1:l
                HbO_sub = HbO((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
                HbR_sub = HbR((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);

                Conc_HbO = Conc_HbO + HbO_sub.*beta_HbO(h);
                Conc_HbR = Conc_HbR + HbR_sub.*beta_HbR(h);
                h = h + 1;
            end
        end
        i_g = i_g + n_g_batch;
    end
    
    Conc.intensity_HbO{i_cond} = Conc_HbO;
    Conc.intensity_HbR{i_cond} = Conc_HbR;
    %% scalp
    
    if cfg.with_scalp == 1
        T_HbO_scalp = T.T_HbO_scalp;
        T_HbR_scalp = T.T_HbR_scalp;
        
        n_g_scalp  = size(G.G_scalp,1);
        tbasis_num = size(T.t_HbO_brain,2);
        n_t_scalp  = tbasis_num*cond_num;
        
        index = 1:n_g_scalp*n_t_scalp;
        index = reshape(index, cond_num*tbasis_num, n_g_scalp);
        index = index((i_cond-1)*tbasis_num+1:i_cond*tbasis_num,:);
        index = reshape(index, [],1);
        
        beta_scalp_HbO = beta(2*n_g*n_t + index); % H_scalp [g1[c1[t1-tn] c2],g2[c1[t1-tn] c2],...]
        beta_scalp_HbR = beta(2*n_g*n_t + n_g_scalp*n_t_scalp + index);

        Conc_HbO = zeros(size(G.G_scalp,2),size(T_HbO_scalp,1),'gpuArray');
        Conc_HbR = zeros(size(G.G_scalp,2),size(T_HbR_scalp,1),'gpuArray');

        n_g_batch = floor(n_g_scalp/n_batch);

        i_g = 1;
        while i_g <= n_g_scalp
            g_start     =   i_g;
            g_end       =   i_g + n_g_batch - 1 ;
            if g_end > n_g_scalp
                g_end = n_g_scalp;
            end
            g_index     =   g_start:g_end;

            Gm = gpuArray(G.G_scalp(g_index,:));
            T_HbO = gpuArray(T_HbO_scalp);
            T_HbR = gpuArray(T_HbR_scalp);
            beta_HbO = beta_scalp_HbO((g_start-1)*tbasis_num+1:g_end*tbasis_num);
            beta_HbR = beta_scalp_HbR((g_start-1)*tbasis_num+1:g_end*tbasis_num);

            HbO = kron(reshape(Gm',[],1),reshape(T_HbO,1,[]));
            HbR = kron(reshape(Gm',[],1),reshape(T_HbR,1,[]));
            [m, n] = size(Gm);
            [k, l] = size(T_HbO);
            h = 1;
            for j = 1:m
                for i = 1:l
                    HbO_sub = HbO((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
                    HbR_sub = HbR((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);

                    Conc_HbO = Conc_HbO + HbO_sub.*beta_HbO(h);
                    Conc_HbR = Conc_HbR + HbR_sub.*beta_HbR(h);
                    h = h + 1;
                end
            end
            i_g = i_g + n_g_batch;
        end

        Conc.intensity_HbO_scalp{i_cond} = Conc_HbO;
        Conc.intensity_HbR_scalp{i_cond} = Conc_HbR;
    end
end
end