function [dc,tHRF,SD, HbO, HbR] = process_sequentially(snirf_file, mat_file, cfg)

Proc_data = process_in_homer_fromOD(snirf_file,mat_file,cfg);

dc   = Proc_data.dcAvg.GetDataTimeSeries('reshape');
tHRF = Proc_data.dcAvg.GetTime();
probe = Proc_data.probe;
SD.Lambda = probe.GetWls();
SD.MeasList = Proc_data.dod.GetMeasList();
SD.SrcPos = probe.sourcePos2D;
SD.DetPos = probe.detectorPos2D;
SD.MeasListAct =  Proc_data.mlActAuto{1,1};
SD.MeasListActAuto =  Proc_data.mlActAuto{1,1};
SD.MeasListActMan =  Proc_data.mlActAuto{1,1};
dc(find(isnan(dc))) = 0;

% get dod conversion for each cond, if more than one condition
dod = [];
for icond = 1:size(dc,4)
    dod(:,:,icond) = hmrConc2OD( squeeze(dc(:,:,:,icond)), SD, [6 6] );
end
n_channels = size(dc,3);
Sensitivity_Matrix = Get_A_dot(cfg.av_data_path, cfg.rhoSD_ssThresh.av,Proc_data.mlActAuto, 1);
Adot = Sensitivity_Matrix.Adot;
Adot_scalp = Sensitivity_Matrix.Adot_scalp;
E = Sensitivity_Matrix.E;
channels = Sensitivity_Matrix.channels;
dod = dod(:,[channels, n_channels+channels],:);

yavgimg = squeeze(dod(find(abs(dod(:,1)) == max(abs(dod(:,1))),1),:,:));

if cfg.AVbrainscalp == 0
    alpha = 0.01;
    Amatrix = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
        squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)];
    [HbO, HbR] = hmrImageReconConc(yavgimg', [], alpha, Amatrix);
    HbO = bsxfun(@rdivide, HbO, Sensitivity_Matrix.ll');
    HbR = bsxfun(@rdivide, HbR, Sensitivity_Matrix.ll');
else
    alpha = 0.01;
    Amatrix = [[squeeze(Adot(:,:,1)) squeeze(Adot_scalp(:,:,1))]*E(1,1) [squeeze(Adot(:,:,1)) squeeze(Adot_scalp(:,:,1))]*E(1,2);
        [squeeze(Adot(:,:,2)) squeeze(Adot_scalp(:,:,2))]*E(2,1) [squeeze(Adot(:,:,2)) squeeze(Adot_scalp(:,:,2))]*E(2,2)];
    [HbO, HbR] = hmrImageReconConc(yavgimg', [], alpha, Amatrix);
    HbO = bsxfun(@rdivide, HbO, [Sensitivity_Matrix.ll Sensitivity_Matrix.ll_scalp]');
    HbR = bsxfun(@rdivide, HbR, [Sensitivity_Matrix.ll Sensitivity_Matrix.ll_scalp]');
end


HbO = HbO(1:20004,:);
HbR = HbR(1:20004,:);