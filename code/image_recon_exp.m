function image_recon_exp(cfg)

cfg.rhoSD_ssThresh.av = 7; % 7 mm

cfg.n_batch_brain = 1;
cfg.n_batch_scalp = 1;
cfg.with_scalp = 0;

At_file = 'atlasViewer.mat';
atlasViewer = load([cfg.av_data_path, At_file]);
warning off

faces = atlasViewer.fwmodel.mesh.faces;
brain_vertices = atlasViewer.fwmodel.mesh.vertices; 

group_folder = dir(fullfile(cfg.datapath,'Subj*'));
subject_folder =  setdiff({group_folder([group_folder.isdir]).name},{'.','..'});

n_vertex_brain = 20004;
image_folder = 'image';
n_subject = numel(subject_folder);

fprintf('In folder %s\n Detect %d subject folders\n', cfg.datapath, n_subject)

data_folder = 'Conc_data';
mkdir(fullfile(cfg.savepath, data_folder));

for ii = 1:n_subject
    snirf_path = fullfile(cfg.datapath, subject_folder{ii});
    snirffile = dir(fullfile(snirf_path,'*.snirf'));
    snirffilelist = {snirffile(~[snirffile.isdir]).name};

    n_snirf = numel(snirffilelist);
    subject_savepath = fullfile(cfg.savepath,'Conc_data');
    if isfile(fullfile(subject_savepath,[subject_folder{ii},'.mat']))
        fprintf('Already exist and loading: %s\n', subject_savepath)
        load(fullfile(subject_savepath,[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
    else
        fprintf('Do not exist and calculating: %s\n', subject_savepath)
        for jj = 1:n_snirf
            snirf_file_name = snirffilelist{jj};
            fprintf('file is %s\n',snirf_file_name)
            calculated_mat = [];
            snirf_file  =   fullfile(snirf_path,snirf_file_name);
            calculated_mat = Proc_stream_v7(snirf_file, [], cfg, calculated_mat);
            calculated_mat.Conc = approximate_Conc(calculated_mat, cfg);
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
        end
        plot_Hb(subject_folder{ii}, faces, brain_vertices, HbO_S, HbR_S, cfg, image_folder,'off')
        save(fullfile(subject_savepath,[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
    end
    
    if ii == 1
        n_cond = size(HbO_S, 2);
        G_HbO_mat = zeros(n_vertex_brain, n_cond, n_subject);
        G_HbR_mat = zeros(n_vertex_brain, n_cond, n_subject);
    end
    G_HbO_mat(:,:,ii) = HbO_S;
    G_HbR_mat(:,:,ii) = HbR_S;
end

% group average
Conc_mean_HbO = mean(G_HbO_mat,3);
Conc_mean_HbR = mean(G_HbR_mat,3);

plot_Hb('group_mean', faces, brain_vertices, Conc_mean_HbO, Conc_mean_HbR, cfg, image_folder,'off')
save(fullfile(cfg.savepath, data_folder, 'group_mean.mat'), 'Conc_mean_HbO','Conc_mean_HbR')

cfg.savepath_results = cfg.savepath;
n_cond = size(G_HbO_mat,2);
for i = 1:n_cond
    output(i) = evaluate_performance_new(G_HbO_mat(:,i,:), G_HbR_mat(:,i,:), mean(G_HbO_mat(:,i,:),3), mean(G_HbR_mat(:,i,:),3), cfg, faces, brain_vertices,i);
end
save(fullfile(cfg.savepath, data_folder, 'output.mat'), 'output')
end