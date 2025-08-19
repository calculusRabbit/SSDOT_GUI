function output = evaluate_performance_new(G_HbO_mat,G_HbR_mat, HbO_S_zero, HbR_S_zero, cfg, faces, brain_vertices, cond)

A_file = ['fw',filesep, 'Adot.mat'];
A_scalp_file = ['fw',filesep,'Adot_scalp.mat'];

Adot        = load([cfg.av_data_path, A_file]);
Adot_scalp  = load([cfg.av_data_path, A_scalp_file]);

M = Make_mask(cfg.mask_threshold, Adot.Adot, Adot_scalp.Adot_scalp);
index = setdiff(1:size(G_HbO_mat,1), M.mask_brain);
G_HbO_mat(index,:,:) = 0;
G_HbR_mat(index,:,:) = 0;
HbO_S_zero(index,:,:) = 0;
HbR_S_zero(index,:,:) = 0;
if exist('cond','var') % this is for experimental data
    if isfile(fullfile(cfg.ROI_path, 'ROI', ['ROI', num2str(cond), '.mat']))
        fprintf('load ROI\n')
        load(fullfile(cfg.ROI_path, 'ROI', ['ROI', num2str(cond), '.mat']), 'index_ROI', 'index_nonROI')
    else
        fprintf('create ROI\n')
        [index_ROI, index_nonROI] = find_ROI(HbO_S_zero,HbR_S_zero, M, brain_vertices, 0.90, 0.01, 40, cfg);
        mkdir(fullfile(cfg.ROI_path, 'ROI'))
        save(fullfile(cfg.ROI_path, 'ROI', ['ROI', num2str(cond), '.mat']), 'index_ROI', 'index_nonROI')
    end
else % this is for simulation data
    if isfile(fullfile(cfg.ROI_path, 'ROI', ['ROI', '.mat']))
        fprintf('load ROI\n')
        load(fullfile(cfg.ROI_path, 'ROI', ['ROI', '.mat']), 'index_ROI', 'index_nonROI')
%     if 0
    else
        fprintf('create ROI\n')
        [index_ROI, index_nonROI] = find_ROI(HbO_S_zero,HbR_S_zero, M, brain_vertices, 0.90, 0.01, 40, cfg);
        mkdir(fullfile(cfg.ROI_path, 'ROI'))
        save(fullfile(cfg.ROI_path, 'ROI', ['ROI', '.mat']), 'index_ROI', 'index_nonROI')
    end
end
% [index_ROI, index_nonROI] = find_ROI(HbO_S_zero,HbR_S_zero, M, brain_vertices, 0.50, 0.01, 0);
axes_order = [2,1,3];

h_figure = figure('name','brain');
intensity = zeros(size(brain_vertices,1), 1);
h = trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','black','facealpha',0.1,'edgealpha',0, 'visible','on');
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);

hold on
plot3(brain_vertices(index_ROI,axes_order(1)),brain_vertices(index_ROI,axes_order(2)),brain_vertices(index_ROI,axes_order(3)),'r.','linewidth',1,'markersize',13);
plot3(brain_vertices(index_nonROI,axes_order(1)),brain_vertices(index_nonROI,axes_order(2)),brain_vertices(index_nonROI,axes_order(3)),'g.','linewidth',1,'markersize',13);
% left
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
% right
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, -2291.8, 130.0])
% camup([-1.0, 0.0, 0.0])

set(gcf,'position',[10         10          313         242])
axis off
axis image
image_folder = 'image';
mkdir(fullfile(cfg.savepath,image_folder))
if exist('cond','var') % this is for experimental data
    savefig(h_figure, fullfile(cfg.savepath,image_folder,['ROI',num2str(cond),'.fig']))
else
    savefig(h_figure, fullfile(cfg.savepath,image_folder,'ROI.fig'))
end
close(h_figure)

%% t value
% group average
Conc_mean_HbO = mean(G_HbO_mat,3);
Conc_mean_HbR = mean(G_HbR_mat,3);

% visualize the subject results
if exist('cond','var') % this is for experimental data
    plot_Hb(['group_mean', num2str(cond)], faces, brain_vertices, Conc_mean_HbO, Conc_mean_HbR, cfg, image_folder,'off')
else
    plot_Hb('group_mean', faces, brain_vertices, Conc_mean_HbO, Conc_mean_HbR, cfg, image_folder,'off')
end

% figure
% subplot(131)
% for subj = 1:size(G_HbO_mat,3)
%     y = G_HbO_mat(:,:,subj);
%     plot(y)
%     hold on
% end
% subplot(132)
% for subj = 1:size(G_HbO_mat,3)
%     y = G_HbO_mat(index_ROI,:,subj);
%     plot(y)
%     hold on
% end
% subplot(133)
% for subj = 1:size(G_HbO_mat,3)
%     y = mean(G_HbO_mat(index_ROI,:,subj));
%     plot(1, y, 'o')
%     hold on
% end
% t value
% first average the voxels
vox_avg_HbO = mean(G_HbO_mat(index_ROI,:,:),1);
vox_avg_HbR = mean(G_HbR_mat(index_ROI,:,:),1);
% then std across subj
std_subj_HbO = squeeze(std(vox_avg_HbO, [], 3,'omitnan'));
std_subj_HbR = squeeze(std(vox_avg_HbR, [], 3,'omitnan'));

mean_subj_HbO = squeeze(mean(vox_avg_HbO, 3,'omitnan'));
mean_subj_HbR = squeeze(mean(vox_avg_HbR, 3,'omitnan'));

t_value_HbO = mean_subj_HbO/(std_subj_HbO/sqrt(size(G_HbO_mat,3)));
t_value_HbR = mean_subj_HbR/(std_subj_HbR/sqrt(size(G_HbR_mat,3)));

output.t_value_HbO = t_value_HbO;
output.t_value_HbR = t_value_HbR;
output.Contrast_HbO = mean_subj_HbO;
output.Contrast_HbR = mean_subj_HbR;
output.std_HbO = std_subj_HbO;
output.std_HbR = std_subj_HbR;

% t value nonROI
% first average the voxels
vox_avg_HbO = mean(G_HbO_mat(index_nonROI,:,:),1);
vox_avg_HbR = mean(G_HbR_mat(index_nonROI,:,:),1);
% then std across subj
std_subj_HbO = squeeze(std(vox_avg_HbO, [], 3,'omitnan'));
std_subj_HbR = squeeze(std(vox_avg_HbR, [], 3,'omitnan'));

mean_subj_HbO = squeeze(mean(vox_avg_HbO, 3,'omitnan'));
mean_subj_HbR = squeeze(mean(vox_avg_HbR, 3,'omitnan'));

t_value_HbO = mean_subj_HbO/(std_subj_HbO/sqrt(size(G_HbO_mat,3)));
t_value_HbR = mean_subj_HbR/(std_subj_HbR/sqrt(size(G_HbR_mat,3)));

output.t_value_HbO2 = t_value_HbO;
output.t_value_HbR2 = t_value_HbR;
output.Contrast_HbO2 = mean_subj_HbO;
output.Contrast_HbR2 = mean_subj_HbR;
output.std_HbO2 = std_subj_HbO;
output.std_HbR2 = std_subj_HbR;

%% test
% thred = 30:2:45;
% t_value_HbO = zeros(size(thred));
% SBR_HbO = zeros(size(thred));
% for i = 1:length(thred)
%     [index_ROI, index_nonROI] = find_ROI(HbO_S_zero,HbR_S_zero, M, brain_vertices, 0.95, 0.001, thred(i));
%     vox_avg_HbO = mean(G_HbO_mat(index_ROI,:,:),1);
%     std_subj_HbO = squeeze(std(vox_avg_HbO, [], 3,'omitnan'));
%     mean_subj_HbO = squeeze(mean(vox_avg_HbO, 3,'omitnan'));
%     t_value_HbO(i) = mean_subj_HbO/(std_subj_HbO/sqrt(size(G_HbO_mat,3)));
%     SBR_HbO(i) = mean(abs(mean(G_HbO_mat(index_ROI,:,:),1)./std(G_HbO_mat(index_nonROI,:,:),[],1)),3);
% end
% figure
% plot(thred, t_value_HbO)
% figure
% plot(thred, SBR_HbO)
%% SBR
SBR_HbO = mean(abs(mean(G_HbO_mat(index_ROI,:,:),1)./std(G_HbO_mat(index_nonROI,:,:),[],1)),3);
SBR_HbR = mean(abs(mean(G_HbR_mat(index_ROI,:,:),1)./std(G_HbR_mat(index_nonROI,:,:),[],1)),3);
% mean_HbO = mean(G_HbO_mat, 3);
% mean_HbR = mean(G_HbR_mat, 3);
% SBR_HbO = mean(mean_HbO(index_ROI, :),1)/std(mean_HbO(index_nonROI, :), [], 1);
% SBR_HbR = mean(mean_HbR(index_ROI, :),1)/std(mean_HbR(index_nonROI, :), [], 1);
output.SBR_HbO = SBR_HbO;
output.SBR_HbR = SBR_HbR;

%% LI
left_hemi_idx  = find(brain_vertices(:,1) > 130);
right_hemi_idx = find(brain_vertices(:,1) <= 130);
LI_subj_HbO = zeros(size(G_HbO_mat, 3), 1);
LI_subj_HbR = zeros(size(G_HbR_mat, 3), 1);
for i = 1:size(G_HbO_mat,3)
    HbO_value = squeeze(G_HbO_mat(:,:,i));
    HbR_value = - squeeze(G_HbR_mat(:,:,i));
    if exist('cond','var')
        LI_subj_HbO(i) = find_LI(HbO_value,left_hemi_idx, right_hemi_idx, cond);
        LI_subj_HbR(i) = find_LI(HbR_value,left_hemi_idx, right_hemi_idx, cond);
    end
end
HbO_value = squeeze(mean(G_HbO_mat, 3));
HbR_value = - squeeze(mean(G_HbR_mat, 3));
if exist('cond','var')
    LI_mean_HbO = find_LI(HbO_value,left_hemi_idx, right_hemi_idx, cond);
    LI_mean_HbR = find_LI(HbR_value,left_hemi_idx, right_hemi_idx, cond);
    output.LI_subj_HbO = LI_subj_HbO;
    output.LI_subj_HbR = LI_subj_HbR;
    output.LI_mean_HbO = LI_mean_HbO;
    output.LI_mean_HbR = LI_mean_HbR;
end

%% positional error
 
if ~exist('cond','var') % this is for simulation data
    intensity = Conc_mean_HbO(M.mask_brain);
    error = 0;
    sum_x = 0;
    axes_order = [2,1,3];
    FVHM_index = find(intensity > 0.5*max(intensity));
%     FVHM_index = index_ROI;
    
    for i = 1:length(FVHM_index)
        index = FVHM_index(i);
        x = intensity(index)-min(intensity(FVHM_index));
        index_brain = M.mask_brain(index);
        r = sqrt((brain_vertices(index_brain,axes_order(1)) - cfg.center(1)).^2 + ...
            (brain_vertices(index_brain,axes_order(2)) - cfg.center(2)).^2 + ...
            (brain_vertices(index_brain,axes_order(3)) - cfg.center(3)).^2);
        error = error + x*r;
        sum_x = sum_x + x;
    end
    if length(FVHM_index) == 1
        output.positional_error_HbO = r;
    else
        output.positional_error_HbO = error/sum_x;
    end
    output.positional_error_HbR = output.positional_error_HbO;
    
%     % second
%     if cfg.original_ROI == 0
%         max_pos_true_HbO = brain_vertices(find(HbO_S_zero == max(HbO_S_zero),1),:);
%         max_pos_true_HbR = brain_vertices(find(HbR_S_zero == max(HbR_S_zero),1),:);
%     elseif cfg.original_ROI == 1
        max_pos_true_HbO = cfg.center(axes_order);
        max_pos_true_HbR = cfg.center(axes_order);
%     end
    % how many max values
    fprintf('max(HbO) %d\t max(HbR) %d \n', length(find(Conc_mean_HbO == max(Conc_mean_HbO),1)), length(find(Conc_mean_HbR == max(Conc_mean_HbR),1)))
    max_pos_HbO = brain_vertices(find(Conc_mean_HbO == max(Conc_mean_HbO),1),:);
    max_pos_HbR = brain_vertices(find(Conc_mean_HbR == min(Conc_mean_HbR),1),:);
    output.positional_error_HbO_2 = sqrt(sum((max_pos_HbO - max_pos_true_HbO).^2));
    output.positional_error_HbR_2 = sqrt(sum((max_pos_HbR - max_pos_true_HbR).^2));
    
end

end

function LI = find_LI(HbO_value,left_hemi_idx, right_hemi_idx, cond)
V_50 = find(HbO_value > 0.5* max(HbO_value));
if cond == 1
    C = length(intersect(V_50, right_hemi_idx));
    I = length(intersect(V_50, left_hemi_idx));
elseif cond == 2
    C = length(intersect(V_50, left_hemi_idx));
    I = length(intersect(V_50, right_hemi_idx));
end
LI = (C - I)/(C + I);
end