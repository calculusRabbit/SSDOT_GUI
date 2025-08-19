function plot_Hb(file, faces, brain_vertices, HbO, HbR, cfg, image_folder, close_fig)
for i_cond = 1:size(HbO,2)
    id_string = sprintf('%s%d',file,i_cond);
    h = figure('Name',id_string);
    if isfield(cfg,'bi') && cfg.bi == 1
        % HbO
        intensity_HbO = HbO(:,i_cond);
        subplot(221)
        plot_intensity(cfg, faces,brain_vertices, intensity_HbO, 'L')
        subplot(222)
        plot_intensity(cfg, faces,brain_vertices, intensity_HbO, 'R')
        % HbR
        intensity_HbR = HbR(:,i_cond);
        subplot(223)
        plot_intensity(cfg, faces,brain_vertices, intensity_HbR, 'L')
        subplot(224)
        plot_intensity(cfg, faces,brain_vertices, intensity_HbR, 'R')
        set(gcf,'position',[ 10         10        895         472])

    else
        % HbO
        intensity_HbO = HbO(:,i_cond);
        subplot(211)
        plot_intensity(cfg, faces,brain_vertices, intensity_HbO, 'L')
        % HbR
        intensity_HbR = HbR(:,i_cond);
        subplot(212)
        plot_intensity(cfg,faces,brain_vertices, intensity_HbR, 'L')
        set(gcf,'position',[ 10         10         399         600])

    end
    mkdir([cfg.savepath,filesep,image_folder])
    savefig(h, fullfile(cfg.savepath,image_folder,[id_string,'.fig']))
    if strcmp(close_fig, 'off')
        close(h)
    end
end
end