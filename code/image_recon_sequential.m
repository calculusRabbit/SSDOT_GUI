function output  = image_recon_sequential(cfg)
cfg.rhoSD_ssThresh.av = 5; % 7 mm
At_file = 'atlasViewer.mat';
atlasViewer = load([cfg.av_data_path, At_file]);
warning off

faces = atlasViewer.fwmodel.mesh.faces;
brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
group_folder = dir(fullfile(cfg.savepath_snirf,'Subj*'));
subject_folder =  setdiff({group_folder([group_folder.isdir]).name},{'.','..'});
n_vertex_brain = 20004;

n_subject = numel(subject_folder);
if strcmp(cfg.savepath(end-4:end), 'brain')
    probeProbePath = fullfile([cfg.savepath,'scalp'],'PlotProbe');
else
    probeProbePath = fullfile(cfg.savepath,'PlotProbe');
end
for ii = 1:n_subject
    if isfile(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},'.mat']))...
            && isfile(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},'_zeros.mat']))...
            && isfile(fullfile(probeProbePath, 'plotProbe.mat'))...
            && isfile(fullfile(probeProbePath,[subject_folder{ii},'.mat']))
        fprintf('load\n')
        load(fullfile(cfg.savepath, 'Conc_data_seq',[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
        load(fullfile(cfg.savepath, 'Conc_data_seq',[subject_folder{ii},'_zeros.mat']), 'HbO_S_zero','HbR_S_zero')
        load(fullfile(probeProbePath, [subject_folder{ii},'.mat']), 'dc_subj')
        if ii == 1
            load(fullfile(probeProbePath,'plotProbe.mat') ,'tHRF','SD')
        end
%     if 0
    else
        snirffile = dir(fullfile(cfg.savepath_snirf,subject_folder{ii},'*.snirf'));
        matfile = dir(fullfile(cfg.savepath_rs,subject_folder{ii},'*.mat'));
        snirffilelist = {snirffile(~[snirffile.isdir]).name};
        matffilelist = {matfile(~[matfile.isdir]).name};
        n_snirf = numel(matffilelist);
        for jj = 1:n_snirf
            file = snirffilelist{1};
            file_mat = matffilelist{jj};
            fprintf('file is %s\n',file)
            if strcmp(cfg.savepath(end-14:end-11), 'noSS') || strcmp(cfg.savepath(end-9:end-6), 'noSS')
                savepath_zero_Hb = cfg.savepath;
            else
                savepath_zero_Hb = cfg.savepath;
                index = strfind(savepath_zero_Hb, filesep);
                index_SS = strfind(savepath_zero_Hb,'SS');

                if strcmp(savepath_zero_Hb(index(3) + 19:index(3) + 23),'gamma')
                    savepath_zero_Hb = [savepath_zero_Hb(1:index(3) + 24),'noSS',savepath_zero_Hb(index_SS + 2:end)];
                else
                    savepath_zero_Hb = [savepath_zero_Hb(1:index(3) + 27),'noSS',savepath_zero_Hb(index_SS + 2:end)];
                end
            end
%             if isfile(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file_mat(1:end-4),'.mat']))...
%                     && isfile(fullfile(savepath_zero_Hb,'Conc_data_seq',['Subject_', num2str(jj), '_zeros.mat']))...
%                     && isfile(fullfile(probeProbePath,'plotProbe.mat') )...
%                     && isfile(fullfile(probeProbePath,[subject_folder{ii},file_mat(1:end-4),'_dc.mat']))
%                 
%                 load(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file_mat(1:end-4),'.mat']), 'HbO','HbR')
%                 load(fullfile(savepath_zero_Hb,'Conc_data_seq',['Subject_', num2str(jj), '_zeros.mat']), 'HbO_zero','HbR_zero')
%                 load(fullfile(probeProbePath,[subject_folder{ii},file_mat(1:end-4),'_dc.mat']),'dc')
%                 if ii == 1 && jj == 1
%                     load(fullfile(probeProbePath,'plotProbe.mat') ,'tHRF','SD')
%                 end
            if 0
            else
                snirf_file = fullfile(cfg.savepath_snirf,subject_folder{ii},file);
                mat_file = fullfile(cfg.savepath_rs,subject_folder{ii},file_mat);
                mat_file_zero = fullfile(cfg.savepath_zeros,'Subject',['Subject_', num2str(jj), '.mat']);
                [dc, tHRF, SD, HbO, HbR] = process_sequentially(snirf_file, mat_file, cfg);
                if strcmp(cfg.savepath(end-14:end-11), 'noSS') || strcmp(cfg.savepath(end-9:end-6), 'noSS')
                    cfg.prune = 0;
                    [~, ~, ~, HbO_zero, HbR_zero] = process_sequentially(snirf_file, mat_file_zero, cfg);
                    cfg.prune = 5;
                    mkdir([cfg.savepath,filesep,'Conc_data_seq'])
                    save(fullfile(savepath_zero_Hb,'Conc_data_seq', ['Subject_', num2str(jj), '_zeros.mat']), 'HbO_zero','HbR_zero')
                else
                    load(fullfile(savepath_zero_Hb,'Conc_data_seq', ['Subject_', num2str(jj), '_zeros.mat']), 'HbO_zero','HbR_zero')
                end
                
                if strcmp(cfg.savepath(end-4:end), 'scalp')
                    % plot probe run level
                    if ii == 1 && jj == 1
                        mkdir(fullfile(cfg.savepath,'PlotProbe'))
                        save(fullfile(cfg.savepath,'PlotProbe','plotProbe.mat'),'tHRF','SD')
                    end
                    save(fullfile(probeProbePath,[subject_folder{ii},file_mat(1:end-4),'_dc.mat']),'dc')
                end
                mkdir([cfg.savepath,filesep,'Conc_data_seq'])
                save(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file_mat(1:end-4),'.mat']), 'HbO','HbR')
            end
            if jj == 1
                HbO_S = HbO./n_snirf;
                HbR_S = HbR./n_snirf;
            else
                HbO_S = HbO_S + HbO./n_snirf;
                HbR_S = HbR_S + HbR./n_snirf;
            end
            if jj == 1
                HbO_S_zero = HbO_zero./n_snirf;
                HbR_S_zero = HbR_zero./n_snirf;
            else
                HbO_S_zero = HbO_S_zero + HbO_zero./n_snirf;
                HbR_S_zero = HbR_S_zero + HbR_zero./n_snirf;
            end
            if jj == 1
                dc_subj = dc./n_snirf;
            else
                dc_subj = dc_subj + dc./n_snirf;
            end
        end
        if strcmp(cfg.savepath(end-4:end), 'scalp')
            save(fullfile(cfg.savepath,'PlotProbe',[subject_folder{ii},'.mat']), 'dc_subj')
        end
        plot_Hb(subject_folder{ii}, faces, brain_vertices, HbO_S, HbR_S, cfg, 'image_seq', 'off')
        save(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
        save(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},'_zeros.mat']), 'HbO_S_zero','HbR_S_zero')
    end
    if ii == 1
        dc_group = dc_subj./n_subject;
    else
        dc_group = dc_group + dc_subj./n_subject;
    end
    
    if ii == 1
        n_cond = size(HbO_S, 2);
        G_HbO_mat = zeros(n_vertex_brain, n_cond, n_subject);
        G_HbR_mat = zeros(n_vertex_brain, n_cond, n_subject);
    end
    G_HbO_mat(:,:,ii) = HbO_S;
    G_HbR_mat(:,:,ii) = HbR_S;
end
if strcmp(cfg.savepath(end-4:end), 'scalp')
    if ~isfile(fullfile(cfg.savepath,'PlotProbe','group.fig'))
        % plotProbe for group
        h = figure;
        plotProbe( dc_group, tHRF, SD, SD, [], [], [], 0, 0);
        axis off
        mkdir(fullfile(cfg.savepath,'PlotProbe'))
        savefig(h, fullfile(cfg.savepath,'PlotProbe','group.fig'))
        close(h)
    end
end
%% scale the image
scaler = max(mean(G_HbO_mat,3));
G_HbO_mat_scaled = G_HbO_mat./scaler;
G_HbR_mat_scaled = G_HbR_mat./scaler;
scaler = max(mean(HbO_S_zero,3));
HbO_S_zero_scaled = HbO_S_zero./scaler;
HbR_S_zero_scaled = HbR_S_zero./scaler;

output = evaluate_performance_new(G_HbO_mat_scaled,G_HbR_mat_scaled, HbO_S_zero_scaled, HbR_S_zero_scaled, cfg, faces, brain_vertices);
save(fullfile(cfg.savepath, 'Conc_data_seq','cfg.mat'), 'cfg', 'output')

end

