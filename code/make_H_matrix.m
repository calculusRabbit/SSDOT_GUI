function [H_brain, H_scalp] = make_H_matrix(A, G, T, Y, device, cfg)

canUseGPU=true;

if strcmp(device, 'gpu')
    try
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
    %% calculate H_brain
    n_batch = cfg.n_batch_brain;
    H_finished = false;
    i = 1;
    while i <= 100
        try
            H_brain = calculate_H_batch_gpu(A.Amatrix_brain, G.G_brain, T.T_HbO_brain, T.T_HbR_brain, Y, n_batch);
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
    %% calculate H_scalp

    n_batch = cfg.n_batch_scalp;
    H_finished = false;
    i = 1;
    while i <= 100
        try
            H_scalp = calculate_H_batch_gpu(A.Amatrix_scalp, G.G_scalp, T.T_HbO_scalp, T.T_HbR_scalp, Y, n_batch);
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

    
end

if strcmp(device, 'cpu') || ~canUseGPU
    [H_brain,H_scalp] = calculate_H_batch_cpu(A, G, T, Y);
end
end

% gpu

function [H, H_ll] = calculate_H_batch_gpu(A, G, T_HbO, T_HbR, Y, n_batch)

n_g = size(G,1);
n_t = size(T_HbO,2);

n_g_batch = floor(n_g/n_batch);
if n_g_batch == 0
    n_g_batch = 1;
end
H = zeros(length(Y),n_g*n_t);

w1_HbO = gpuArray(A.w1_HbO);
w1_HbR = gpuArray(A.w1_HbR);
w2_HbO = gpuArray(A.w2_HbO);
w2_HbR = gpuArray(A.w2_HbR);

i_g = 1;
while i_g <= n_g
    g_start     =   i_g;
    g_end       =   i_g + n_g_batch - 1 ;
    if g_end > n_g
        g_end = n_g;
    end
    g_index     =   g_start:g_end;
    
    G_brain = gpuArray(G(g_index,:));
    T_HbO = gpuArray(T_HbO);
    T_HbR = gpuArray(T_HbR);
    
    HbO  = kron(reshape(G_brain',[],1),reshape(T_HbO,1,[]));
    HbR  = kron(reshape(G_brain',[],1),reshape(T_HbR,1,[]));
    [m, n] = size(G_brain);
    [k, l] = size(T_HbO);
    for j = 1:m
        for i = 1:l
            HbO_sub = HbO((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            HbR_sub = HbR((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            Agt_w1_HbO = w1_HbO*HbO_sub;
            Agt_w1_HbR = w1_HbR*HbR_sub;
            Agt_w2_HbO = w2_HbO*HbO_sub;
            Agt_w2_HbR = w2_HbR*HbR_sub;
            
            Agt_w1_HbO = gather(Agt_w1_HbO);
            Agt_w1_HbR = gather(Agt_w1_HbR);
            Agt_w2_HbO = gather(Agt_w2_HbO);
            Agt_w2_HbR = gather(Agt_w2_HbR);
            
            H(1:length(Y)/2,(g_start + j- 2) * l + i ) = reshape(Agt_w1_HbO,[],1);
            H(1:length(Y)/2, n_g*n_t + (g_start + j- 2) * l + i ) = reshape(Agt_w1_HbR,[],1);
            H(length(Y)/2+1:end,(g_start + j- 2) * l + i ) = reshape(Agt_w2_HbO,[],1);
            H(length(Y)/2+1:end,n_g*n_t + (g_start + j- 2) * l + i ) = reshape(Agt_w2_HbR,[],1);
        end
        
    end
    
    i_g = i_g + n_g_batch;
end

end

% cpu
function [H_brain,H_scalp] = calculate_H_batch_cpu(A, G, T, Y)

%% H_brain
n_g = size(G.G_brain,1);
n_t = size(T.T_HbO_brain,2);
H_brain = zeros(length(Y),n_g*n_t*2);
j = 1;
for i_g = 1:n_g
    for i_t = 1:n_t
        gt_HbO = kron(G.G_brain(i_g,:),T.T_HbO_brain(:,i_t));
        gt_HbR = kron(G.G_brain(i_g,:),T.T_HbR_brain(:,i_t));
        Agt_w1_HbO = A.Amatrix_brain.w1_HbO*gt_HbO';
        Agt_w1_HbR = A.Amatrix_brain.w1_HbR*gt_HbR';
        Agt_w2_HbO = A.Amatrix_brain.w2_HbO*gt_HbO';
        Agt_w2_HbR = A.Amatrix_brain.w2_HbR*gt_HbR'; % channels X time 
        
        H_brain(1:length(Y)/2,j)                = reshape(Agt_w1_HbO,[],1);
        H_brain(1:length(Y)/2, n_g*n_t + j)     = reshape(Agt_w1_HbR,[],1);
        H_brain(length(Y)/2+1:end,j)            = reshape(Agt_w2_HbO,[],1);
        H_brain(length(Y)/2+1:end, n_g*n_t + j) = reshape(Agt_w2_HbR,[],1);
        %         HbO                |     HbR
        % w1 | c1[t1-tn] c2 c1 c2 ...       | c1 c2 c1 c2...  
        % w2 | c1[t1-tn] c2 c1 c2 ...       | c1 c2 c1 c2...
        j = j + 1;
    end
end
%% H_scalp
n_g = size(G.G_scalp,1);
n_t = size(T.T_HbO_scalp,2);
H_scalp = zeros(length(Y),n_g*n_t*2);

j = 1;
for i_g = 1:n_g
    for i_t = 1:n_t
        gt_HbO = kron(G.G_scalp(i_g,:),T.T_HbO_scalp(:,i_t));
        gt_HbR = kron(G.G_scalp(i_g,:),T.T_HbR_scalp(:,i_t));
        Agt_w1_HbO = A.Amatrix_scalp.w1_HbO*gt_HbO';
        Agt_w1_HbR = A.Amatrix_scalp.w1_HbR*gt_HbR';
        Agt_w2_HbO = A.Amatrix_scalp.w2_HbO*gt_HbO';
        Agt_w2_HbR = A.Amatrix_scalp.w2_HbR*gt_HbR';
        
        H_scalp(1:length(Y)/2,j)                = reshape(Agt_w1_HbO,[],1);
        H_scalp(1:length(Y)/2, n_g*n_t + j)     = reshape(Agt_w1_HbR,[],1);
        H_scalp(length(Y)/2+1:end,j)            = reshape(Agt_w2_HbO,[],1);
        H_scalp(length(Y)/2+1:end, n_g*n_t + j) = reshape(Agt_w2_HbR,[],1);
        
        j = j + 1;
    end
end


end