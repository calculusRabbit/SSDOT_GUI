function output = reconstruct_simulate_L_and_CNR_rs_load(cfg)
n = length(cfg.alpha_list);
x_mat = zeros(n, 1);
Ax_y_mat = zeros(n, 1);
Contrast = zeros(n, 1);
image_CNR = zeros(n, 1);
error = zeros(n, 1);

for i = 1:n
    alpha = cfg.alpha_list(i);
    fprintf('%d, alpha = %d\n',i, alpha)
    rs_folder = fullfile(cfg.datapath,['results_',sprintf('%.3e',alpha)],'Conc_data');
    load(fullfile(rs_folder,'cfg.mat'), 'output', 'abs_x_G', 'Ax_y_G', 'CNR_G')
    x_mat(i) = abs_x_G;
    Ax_y_mat(i) = Ax_y_G;
    Contrast(i) = CNR_G.Contrast;
    image_CNR(i) = CNR_G.image_CNR;
    error(i) = CNR_G.error;
end

curve.x_mean = x_mat;
curve.Ax_y_mean = Ax_y_mat;
curve.Contrast_mean = Contrast;
curve.image_CNR_mean = image_CNR;
curve.error_mean = error;

plot_curves(curve,cfg)
end


