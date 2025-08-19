function plot_simulation(cfg)
light_red = [247 192 192]./255;
rs_folder = dir(fullfile(cfg.savepath,'simulated_rs', num2str(cfg.scale_factor), 'Subject*'));
subject_folder = setdiff({rs_folder([rs_folder.isdir]).name},{'.','..'});
for ii = 4:4
    rs_mat_path = dir(fullfile(cfg.savepath,'simulated_rs', num2str(cfg.scale_factor), subject_folder{ii}, '*.mat'));
    rs_mat_filelist = {rs_mat_path(~[rs_mat_path.isdir]).name};
    %% resting data
    resting_mat_file = fullfile(cfg.savepath, 'resting', [subject_folder{ii},'.mat']);
    load(resting_mat_file,'OD_rearrange','stim','t_new')
    t = t_new;
    
    h = figure(); 
    tlo = tiledlayout(h,4,1);
    nexttile(tlo)

    plot(t, OD_rearrange(:,83)', 'color', light_red)
    xlim([min(t),max(t)])
    ylim([min(OD_rearrange(:,83)), max(OD_rearrange(:, 83))])
    yticks([-0.05 0.05])
    box off
    set(gca,'xtick',[])
    set(gca,'fontname','Arial')
    set(gca, 'fontsize',12)
    
    y = OD_rearrange(:, 83);
    lpf = 0.5;
    fs = 1/(t(2) - t(1));
    lpf_norm = lpf / (fs / 2);
    FilterOrder = 3;
    [z, p, k] = butter(FilterOrder, lpf_norm, 'low');
    [sos, g] = zp2sos(z, p, k);
    OD_rearrange_filtered = filtfilt(sos, g, double(y)); 
    
    for i = 1:length(rs_mat_filelist)
        load(fullfile(cfg.savepath,'simulated_rs', num2str(cfg.scale_factor), subject_folder{ii}, rs_mat_filelist{i}), 'OD');
        
        s_array = stim{i};
        s_index = find(s_array==1);
        y = OD(:, 84);
        lpf = 0.5;
        fs = 1/(t(2) - t(1));
        lpf_norm = lpf / (fs / 2);
        FilterOrder = 3;
        [z, p, k] = butter(FilterOrder, lpf_norm, 'low');
        [sos, g] = zp2sos(z, p, k);
        OD_filtered = filtfilt(sos, g, double(y)); 
         
        nexttile(tlo)
        hold on;
        n_s = length(s_index);
        arrayfun(@(i)xline(t(s_index(i)),':', 'color', [120 120 120]./255, 'linewidth',1.5), 1:n_s)
        
        plot(t, OD_rearrange(:, 83), 'color', light_red);
        plot(t, OD_rearrange_filtered, 'color', 'r', 'linewidth',1.5);
        plot(t, OD_filtered, 'color', [0 0 255]./255, 'linewidth',1.5);
        
        yticks([-0.05 0.05])
        hold off;
        ylim([min([OD_rearrange(:, 83);OD_rearrange_filtered; OD_filtered]), max([OD_rearrange(:, 83); OD_rearrange_filtered; OD_filtered])])
        xlim([min(t),max(t)])
        set(gca, 'fontsize',12)
        box off
        set(gca,'xtick',[100 200 300])
        set(gca,'fontname','Arial')
    end
    ylabel(tlo, '\Delta OD')
    xlabel(tlo, 'Seconds')
    
    set(gcf,'position', [20         20        841         474])
    mkdir(fullfile(cfg.savepath,'simulated_results','plot_simulation', num2str(cfg.scale_factor)));
    savefig(h, fullfile(cfg.savepath,'simulated_results','plot_simulation', num2str(cfg.scale_factor),['plot_simlation_', subject_folder{ii}, '.fig']))
end