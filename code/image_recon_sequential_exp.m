function image_recon_sequential_exp(cfg)

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
            && isfile(fullfile(probeProbePath, 'plotProbe.mat'))...
            && isfile(fullfile(probeProbePath,[subject_folder{ii},'.mat']))
        
        load(fullfile(cfg.savepath, 'Conc_data_seq',[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
        load(fullfile(probeProbePath, [subject_folder{ii},'.mat']), 'dc_subj')
        if ii == 1
            load(fullfile(probeProbePath,'plotProbe.mat') ,'tHRF','SD')
        end
    else
        snirffile = dir(fullfile(cfg.savepath_snirf,subject_folder{ii},'*.snirf'));
        snirffilelist = {snirffile(~[snirffile.isdir]).name};
        n_snirf = numel(snirffilelist);
        for jj = 1:n_snirf
            file = snirffilelist{jj};
            fprintf('file is %s\n',file)
            if isfile(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file(1:end-6),'.mat']))...
                    && isfile(fullfile(probeProbePath,'plotProbe.mat'))...
                    && isfile(fullfile(probeProbePath,[subject_folder{ii},file(1:end-6),'_dc.mat']))
                
                load(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file(1:end-6),'.mat']), 'HbO','HbR')
                load(fullfile(probeProbePath,[subject_folder{ii},file(1:end-6),'_dc.mat']),'dc')
                if ii == 1 && jj == 1
                    load(fullfile(probeProbePath,'plotProbe.mat') ,'tHRF','SD')
                end
            else
                snirf_file = fullfile(cfg.savepath_snirf,subject_folder{ii},file);
                [dc, tHRF, SD, HbO, HbR] = process_sequentially(snirf_file, [], cfg);
                
                if strcmp(cfg.savepath(end-4:end), 'scalp')
                    % plot probe run level
                    if ii == 1 && jj == 1
                        mkdir(fullfile(cfg.savepath,'PlotProbe'))
                        save(fullfile(cfg.savepath,'PlotProbe','plotProbe.mat'),'tHRF','SD')
                    end
                end
                mkdir([cfg.savepath,filesep,'Conc_data_seq'])
                save(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},file(1:end-6),'.mat']), 'HbO','HbR')
                
            end
            if jj == 1
                HbO_S = HbO./n_snirf;
                HbR_S = HbR./n_snirf;
            else
                HbO_S = HbO_S + HbO./n_snirf;
                HbR_S = HbR_S + HbR./n_snirf;
            end
            
            if jj == 1
                dc_subj = dc./n_snirf;
            else
                dc_subj = dc_subj + dc./n_snirf;
            end
        end
        if strcmp(cfg.savepath(end-4:end), 'scalp')
            % plotProbe for this subject
            for cond = 1:size(dc_subj, 4)
                h = figure;
                plotProbe( dc_subj(:,:,:,cond), tHRF, SD, SD, [], [0.5 1.5], [], 0, 0);
                axis off
                set(gcf, 'position', [10         10         397         430])
                savefig(h, fullfile(cfg.savepath,'PlotProbe',[subject_folder{ii},int2str(cond), '.fig']))
                close(h)
            end
            save(fullfile(cfg.savepath,'PlotProbe',[subject_folder{ii}, '.mat']), 'dc_subj')
        end
        % image recon for subject level
        plot_Hb(subject_folder{ii}, faces, brain_vertices, HbO_S, HbR_S, cfg, 'image_seq','off')
        save(fullfile(cfg.savepath,'Conc_data_seq',[subject_folder{ii},'.mat']), 'HbO_S','HbR_S')
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
    for cond = 1:size(dc_group, 4)
        if ~isfile(fullfile(cfg.savepath,'PlotProbe',['group_',int2str(cond), '.fig']))
            h = figure;
            plotProbe( dc_group(:,:,:,cond), tHRF, SD, SD, [], [0.5 1.5], [], 0, 0);
            axis off
            set(gcf, 'position', [10         10         397         430])
            mkdir(fullfile(cfg.savepath,'PlotProbe'))
            savefig(h, fullfile(cfg.savepath,'PlotProbe',['group_mean',int2str(cond), '.fig']))
            close(h)
        end
    end

end
% group average
Conc_mean_HbO = mean(G_HbO_mat,3);
Conc_mean_HbR = mean(G_HbR_mat,3);

plot_Hb('group_mean', faces, brain_vertices, Conc_mean_HbO, Conc_mean_HbR, cfg, 'image_seq', 'off')
save(fullfile(cfg.savepath, 'Conc_data_seq', 'group_mean.mat'), 'Conc_mean_HbO','Conc_mean_HbR')

cfg.savepath_results = cfg.savepath;
n_cond = size(G_HbO_mat,2);
for i = 1:n_cond
    output(i) = evaluate_performance_new(G_HbO_mat(:,i,:), G_HbR_mat(:,i,:), mean(G_HbO_mat(:,i,:),3), mean(G_HbR_mat(:,i,:),3), cfg, faces, brain_vertices,i);
end
save(fullfile(cfg.savepath, 'Conc_data_seq', 'output.mat'), 'output')
end

