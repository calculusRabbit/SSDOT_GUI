function image_recon(cfg)

cfg.rhoSD_ssThresh.av = 7; % 7 mm

cfg.n_batch_brain = 1;
cfg.n_batch_scalp = 1;
cfg.with_scalp = 0;

% acquired_path.av = cfg.av_data_path;
At_file = 'atlasViewer.mat';
atlasViewer = load([cfg.av_data_path, At_file]);
warning off

faces = atlasViewer.fwmodel.mesh.faces;
brain_vertices = atlasViewer.fwmodel.mesh.vertices; 
group_folder = dir(fullfile(cfg.savepath_rs,'Subj*'));
subject_folder =  setdiff({group_folder([group_folder.isdir]).name},{'.','..'});

n_vertex_brain = 20004;
image_folder = 'image';
n_subject = numel(subject_folder);
fprintf('In folder %s\n Detect %d subject folders\n', cfg.savepath_rs, n_subject)

abs_x_G = 0;
Ax_y_G = 0;
CNR_G.Contrast = 0;
CNR_G.image_CNR = 0;
CNR_G.error = 0;
mkdir([cfg.savepath, filesep, 'Conc_data'])

for ii = 1:n_subject
    snirf_path = fullfile(cfg.savepath_snirf, subject_folder{ii});
    snirffile = dir(fullfile(snirf_path,'*.snirf'));
    snirffilelist = {snirffile(~[snirffile.isdir]).name};
    
    mat_path = fullfile(cfg.savepath_rs, subject_folder{ii});
    matfile = dir(fullfile(mat_path,'*.mat'));
    matfilelist = {matfile(~[matfile.isdir]).name};
    
    mat_path_zero = fullfile(cfg.savepath_zeros, 'Subject');
    
    n_snirf = numel(matfilelist);
    subject_savepath = fullfile(cfg.savepath,'Conc_data');
    if isfile(fullfile(subject_savepath,[subject_folder{ii},'.mat'])) ...
            && isfile(fullfile(subject_savepath,'group_zero.mat'))...
            && isfile(fullfile(subject_savepath,[subject_folder{ii},'_cnr.mat']))
        
        fprintf('Already exist and loading: %s\n', subject_savepath)
        load(fullfile(subject_savepath,[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
        load(fullfile(subject_savepath,'group_zero.mat'), 'HbO_S_zero', 'HbR_S_zero')
        load(fullfile(subject_savepath,[subject_folder{ii},'_cnr.mat']), 'abs_x_S', 'Ax_y_S', 'CNR_S')
    else
        fprintf('Do not exist and calculating: %s\n', subject_savepath)
        abs_x_S = 0;
        Ax_y_S = 0;
        CNR_S.Contrast = 0;
        CNR_S.image_CNR = 0;
        CNR_S.error = 0;
        % zeros
        SS_flag = cfg.SS_flag;
        drift_flag = cfg.drift_flag;
        prune = cfg.prune;
%         spatially_regu = cfg.spatially_regu;
        
%         cfg.spatially_regu = 0;
        cfg.SS_flag = 0;
        cfg.drift_flag = 0;
        cfg.prune = 0;
        if isfile(fullfile(subject_savepath,'group_zero.mat'))
            load(fullfile(subject_savepath,'group_zero.mat'), 'HbO_S_zero', 'HbR_S_zero')
        else
            snirf_file_name = snirffilelist{1};
            for jj = 1:n_snirf
                mat_file_name_zero = ['Subject_', num2str(jj), '.mat'];
                fprintf('zero file is %s\n',mat_file_name_zero)
                calculated_mat_zero = [];
                snirf_file  =   fullfile(snirf_path,snirf_file_name);
                mat_file_zero    =   fullfile(mat_path_zero,  mat_file_name_zero);
                calculated_mat_zero = Proc_stream_v7(snirf_file, mat_file_zero, cfg, calculated_mat_zero);
                calculated_mat_zero.Conc = approximate_Conc(calculated_mat_zero, cfg);
                [HbO_zero, HbR_zero] = project2vertex(calculated_mat_zero, n_vertex_brain);
                if jj == 1
                    HbO_S_zero = HbO_zero./n_snirf;
                    HbR_S_zero = HbR_zero./n_snirf;
                else
                    HbO_S_zero = HbO_S_zero + HbO_zero./n_snirf;
                    HbR_S_zero = HbR_S_zero + HbR_zero./n_snirf;
                end
            end
            if cfg.spatially_regu
                for cond = 1:size(HbO_S_zero, 2)
                    HbO_S_zero(:,cond) = bsxfun(@rdivide, HbO_S_zero(:,cond), calculated_mat_zero.ll');
                    HbR_S_zero(:,cond) = bsxfun(@rdivide, HbR_S_zero(:,cond), calculated_mat_zero.ll');
                end
            end
            
            plot_Hb('group_zero', faces, brain_vertices, HbO_S_zero, HbR_S_zero, cfg, image_folder,'off')
            save(fullfile(subject_savepath,'group_zero.mat'), 'HbO_S_zero','HbR_S_zero')
        end
        cfg.SS_flag = SS_flag;
        cfg.drift_flag = drift_flag;
        cfg.prune = prune;
%         cfg.spatially_regu = spatially_regu;
        
        snirf_file_name = snirffilelist{1};
        for jj = 1:n_snirf
            mat_file_name = matfilelist{jj};
            fprintf('file is %s\n',mat_file_name)
            calculated_mat = [];

            snirf_file  =   fullfile(snirf_path,snirf_file_name);
            mat_file    =   fullfile(mat_path,  mat_file_name);
            
            calculated_mat = Proc_stream_v7(snirf_file, mat_file, cfg, calculated_mat);
            calculated_mat.Conc = approximate_Conc(calculated_mat, cfg);
            [abs_x,Ax_y] = calculate_x(calculated_mat);
            CNR = calculate_cnr(calculated_mat, cfg, calculated_mat.Conc);
            [HbO, HbR] = project2vertex(calculated_mat, n_vertex_brain);
            
            if jj == 1
                HbO_S = HbO./n_snirf;
                HbR_S = HbR./n_snirf;
            else
                HbO_S = HbO_S + HbO./n_snirf;
                HbR_S = HbR_S + HbR./n_snirf;
            end
            if cfg.spatially_regu
                for cond = 1:size(HbO, 2)
                    HbO_S(:,cond) = bsxfun(@rdivide, HbO_S(:,cond), calculated_mat.ll');
                    HbR_S(:,cond) = bsxfun(@rdivide, HbR_S(:,cond), calculated_mat.ll');
                end
            end
            abs_x_S = abs_x_S + abs_x./n_snirf;
            Ax_y_S = Ax_y_S + Ax_y./n_snirf;
            CNR_S.Contrast = CNR_S.Contrast + CNR.Contrast./n_snirf;
            CNR_S.image_CNR = CNR_S.image_CNR + CNR.image_CNR./n_snirf;
            CNR_S.error = CNR_S.error + CNR.error./n_snirf;
        end
        % visualize the subject results
        plot_Hb(subject_folder{ii}, faces, brain_vertices, HbO_S, HbR_S, cfg, image_folder,'off')
        % save the subject level results
        save(fullfile(subject_savepath,[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
        save(fullfile(subject_savepath,[subject_folder{ii},'_cnr.mat']), 'abs_x_S', 'Ax_y_S', 'CNR_S')
    end
    
    if ii == 1
        n_cond = size(HbO_S, 2);
        G_HbO_mat = zeros(n_vertex_brain, n_cond, n_subject);
        G_HbR_mat = zeros(n_vertex_brain, n_cond, n_subject);
    end
    G_HbO_mat(:,:,ii) = HbO_S;
    G_HbR_mat(:,:,ii) = HbR_S;
    
    abs_x_G = abs_x_G + abs_x_S./n_subject;
    Ax_y_G = Ax_y_G + Ax_y_S./n_subject;
    CNR_G.Contrast = CNR_G.Contrast + CNR_S.Contrast./n_subject;
    CNR_G.image_CNR = CNR_G.image_CNR + CNR_S.image_CNR./n_subject;
    CNR_G.error = CNR_G.error + CNR_S.error./n_subject;
end
%% scale the image
scaler = max(mean(G_HbO_mat,3));
G_HbO_mat_scaled = G_HbO_mat./scaler;
G_HbR_mat_scaled = G_HbR_mat./scaler;
scaler = max(mean(HbO_S_zero,3));
HbO_S_zero_scaled = HbO_S_zero./scaler;
HbR_S_zero_scaled = HbR_S_zero./scaler;

output = evaluate_performance_new(G_HbO_mat_scaled, G_HbR_mat_scaled, HbO_S_zero_scaled, HbR_S_zero_scaled, cfg, faces, brain_vertices);
% how i get the centrio
% mean_image = mean(G_HbO_mat_scaled,3);
% find(mean_image == max(mean_image))
% 17602


CNR_G.Contrast = gather(CNR_G.Contrast);
CNR_G.image_CNR = gather(CNR_G.image_CNR);
CNR_G.error = gather(CNR_G.error);
abs_x_G = gather(abs_x_G);
Ax_y_G = gather(Ax_y_G);

save(fullfile(cfg.savepath, 'Conc_data','cfg.mat'), 'cfg', 'output', 'abs_x_G', 'Ax_y_G', 'CNR_G')

end