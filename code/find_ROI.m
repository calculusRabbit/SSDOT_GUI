function [index_ROI, index_nonROI] = find_ROI(HbO_S_zero,HbR_S_zero, M, brain_vertices, max_thred, min_thred, distance_thred, cfg)

if cfg.original_ROI == 0
    intensity = HbO_S_zero;
    index_ROI = find(intensity > max_thred*max(intensity));
    index_ROI = intersect(index_ROI, M.mask_brain);
    index_nonROI = find(abs(intensity) < min_thred*max(intensity));
    for i = length(index_nonROI):-1:1
        for j = 1:length(index_ROI)
            distance = sqrt(sum((brain_vertices(index_nonROI(i),:) - brain_vertices(index_ROI(j),:)).^2));
    %         fprintf('%f\n', distance)
            if distance < distance_thred
                index_nonROI(i) = [];
                break
            end
        end
    end

    index_nonROI = intersect(index_nonROI, M.mask_brain);
else

    axes_order = [2,1,3];
    index_perturbation = [];

    for i = 1:size(cfg.center,1)
        perturbation_center = cfg.center(i,:);
        distance = sqrt((brain_vertices(:,axes_order(1)) - perturbation_center(1)).^2 + ...
            (brain_vertices(:,axes_order(2)) - perturbation_center(2)).^2 + ...
            (brain_vertices(:,axes_order(3)) - perturbation_center(3)).^2);
        index_perturbation = [index_perturbation; find(distance<15)];
    end
    index_ROI = index_perturbation;
    index_nonROI = setdiff(M.mask_brain, index_perturbation);
end

if isempty(index_ROI)
    fprintf('index ROI is empty\n')
end
if isempty(index_nonROI)
    fprintf('index nonROI is empty\n')
end