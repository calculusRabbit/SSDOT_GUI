function [HTY, HTH] = make_HTY_and_HTH(H_brain,H_scalp, OD_SS, OD_drift, Y, chunk, device, R)
% H = [H_brain, H_scalp, OD_SS, OD_drift]
canUseGPU=true;

if strcmp(device,'gpu')
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
    
    H_finished = false;
    k = 1;
    while k <= 100
        try
            [HTY, HTH] = calculate_HTY_and_HTH_batch_gpu(H_brain, H_scalp, OD_SS, OD_drift,  Y, chunk, R);
            H_finished = true;
            break
        catch
            chunk = chunk / 2;
            fprintf('chunk is %d\n',chunk)
            k = k + 1;
            continue
        end
    end
    if H_finished == false
        fprintf('Cannot fit into GPU memory even with chuck = %d \n', chunk)
        return
    end
    
end
if strcmp(device, 'cpu') || ~canUseGPU
    H = [H_brain, H_scalp, OD_SS, OD_drift]; 
    HTY = (H'*Y);
    HTH = H'*H;
end

end

function [HTY, HTH] = calculate_HTY_and_HTH_batch_gpu(H_brain,H_scalp,  OD_SS, OD_drift, Y, chunk, R)
R = gpuArray(R);
Y = gpuArray(Y);
HTH = zeros(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + size(OD_drift,2));
HTY = zeros(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + size(OD_drift,2),1); 
% chunk = 100;
for i = 1:chunk:size(H_brain,2)
    index_start = i;
    index_end = i + chunk - 1;
    if index_end > size(H_brain,2)
        index_end = size(H_brain,2);
    end
    index1 = index_start : index_end;
    H_chunk_1 = gpuArray(H_brain(:,index1));
    if ~isempty(R)
        hty = H_chunk_1' * R * Y;
    else
        hty = H_chunk_1'*Y;
    end
    HTY(index1,:) = hty;
    for j = 1:chunk:size(H_brain,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(H_brain,2)
            index_end = size(H_brain,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(H_brain(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(index1,index2) = hth;
    end
    for j = 1:chunk:size(H_scalp,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(H_scalp,2)
            index_end = size(H_scalp,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(H_scalp(:,index2));
        
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        
        HTH(index1,size(H_brain,2) + index2) = hth;
        HTH(size(H_brain,2) + index2,index1) = hth';
    end
    if ~isempty(OD_SS)
        for j = 1:chunk:size(OD_SS,2)
            index_start = j;
            index_end = j + chunk - 1;
            if index_end > size(OD_SS,2)
                index_end = size(OD_SS,2);
            end
            index2 = index_start : index_end;
            H_chunk_2 = gpuArray(OD_SS(:,index2));
            if ~isempty(R)
                hth = H_chunk_1' * R * H_chunk_2;
            else
                hth = H_chunk_1'*H_chunk_2;
            end
            HTH(index1,size(H_brain,2) + size(H_scalp,2)+ index2) = hth;
            HTH(size(H_brain,2) + size(H_scalp,2) + index2,index1) = hth';
        end
    end
    for j = 1:chunk:size(OD_drift,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(OD_drift,2)
            index_end = size(OD_drift,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(OD_drift(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(index1,size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2)+ index2) = hth;
        HTH(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + index2,index1) = hth';
    end
end
for i = 1:chunk:size(H_scalp,2)
    index_start = i;
    index_end = i + chunk - 1;
    if index_end > size(H_scalp,2)
        index_end = size(H_scalp,2);
    end
    index1 = index_start : index_end;
    H_chunk_1 = gpuArray(H_scalp(:,index1));
    if ~isempty(R)
        hty = H_chunk_1' * R * Y;
    else
        hty = H_chunk_1'*Y;
    end
    HTY(size(H_brain,2) + index1,:) = hty;
    for j = 1:chunk:size(H_scalp,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(H_scalp,2)
            index_end = size(H_scalp,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(H_scalp(:,index2));
        if ~isempty(R)
            hth = H_chunk_1'*R*H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(size(H_brain,2) + index1 , size(H_brain,2) + index2) = hth;
    end
    if ~isempty(OD_SS)
        for j = 1:chunk:size(OD_SS,2)
            index_start = j;
            index_end = j + chunk - 1;
            if index_end > size(OD_SS,2)
                index_end = size(OD_SS,2);
            end
            index2 = index_start : index_end;
            H_chunk_2 = gpuArray(OD_SS(:,index2));
            if ~isempty(R)
                hth = H_chunk_1' * R * H_chunk_2;
            else
                hth = H_chunk_1'*H_chunk_2;
            end
            HTH(size(H_brain,2) + index1,size(H_brain,2) + size(H_scalp,2)+ index2) = hth;
            HTH(size(H_brain,2) + size(H_scalp,2) + index2, size(H_brain,2) + index1) = hth';
        end
    end
    for j = 1:chunk:size(OD_drift,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(OD_drift,2)
            index_end = size(OD_drift,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(OD_drift(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(size(H_brain,2) + index1,size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2)+ index2) = hth;
        HTH(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + index2,size(H_brain,2) + index1) = hth';
    end
end
for i = 1:chunk:size(OD_SS,2)
    index_start = i;
    index_end = i + chunk - 1;
    if index_end > size(OD_SS,2)
        index_end = size(OD_SS,2);
    end
    index1 = index_start : index_end;
    H_chunk_1 = gpuArray(OD_SS(:,index1));
    if ~isempty(R)
        hty = H_chunk_1' * R * Y;
    else
        hty = H_chunk_1' * Y;
    end
    HTY(size(H_brain,2) + size(H_scalp,2) + index1,:) = hty;
    for j = 1:chunk:size(OD_SS,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(OD_SS,2)
            index_end = size(OD_SS,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(OD_SS(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(size(H_brain,2) + size(H_scalp,2) + index1 , size(H_brain,2)+ size(H_scalp,2) + index2) = hth;
    end
    for j = 1:chunk:size(OD_drift,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(OD_drift,2)
            index_end = size(OD_drift,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(OD_drift(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R * H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        HTH(size(H_brain,2) + size(H_scalp,2) + index1, size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + index2) = hth;
        HTH(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + index2, size(H_brain,2) + size(H_scalp,2) + index1) = hth';
    end
    
end
for i = 1:chunk:size(OD_drift,2)
    index_start = i;
    index_end = i + chunk - 1;
    if index_end > size(OD_drift,2)
        index_end = size(OD_drift,2);
    end
    index1 = index_start : index_end;
    H_chunk_1 = gpuArray(OD_drift(:,index1));
    if ~isempty(R)
        hty = H_chunk_1'*R*Y;
    else
        hty = H_chunk_1'*Y;
    end
    HTY(size(H_brain,2) + size(H_scalp,2) + size(OD_SS,2) + index1,:) = hty;
    for j = 1:chunk:size(OD_drift,2)
        index_start = j;
        index_end = j + chunk - 1;
        if index_end > size(OD_drift,2)
            index_end = size(OD_drift,2);
        end
        index2 = index_start : index_end;
        H_chunk_2 = gpuArray(OD_drift(:,index2));
        if ~isempty(R)
            hth = H_chunk_1' * R *H_chunk_2;
        else
            hth = H_chunk_1'*H_chunk_2;
        end
        
        HTH(size(H_brain,2) + size(H_scalp,2)+ size(OD_SS,2) + index1 , size(H_brain,2)+ size(H_scalp,2)+ size(OD_SS,2) + index2) = hth;
    end
end
end
