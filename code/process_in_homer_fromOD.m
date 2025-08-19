function Proc_data = process_in_homer_fromOD(file, mat_file,cfg)

acquired = SnirfClass(file);
data = acquired.data;
probe = acquired.probe;

mlActMan =  cell(1,1);
mlActMan{1,1} = ones(length(data.measurementList),1);

tIncMan =  cell(1,1);
tIncMan{1,1} = ones(size(acquired.data.time));


Aaux = [];
tIncAuto = [];
rcMap = [];
stim = acquired.stim;
%% homer3 GLM processing
if isfield(cfg,'prune')
    prune_parameter = cfg.prune;
else
    prune_parameter = 0;
end
prune_SD = cfg.prune_SD;
% prune_parameter = 5;
%% delete other conditions
if isfield(cfg,'condition')
    index = [];
    for i = 1:length(stim)
        stim_i = stim(i);
        if ~strcmp(stim_i.name,cfg.condition)
            index = [index;i];
        end
    end
    stim(index) = [];
end

mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,cfg.prune_range,prune_parameter,prune_SD);
dod = hmrR_Intensity2OD(data);
if ~isempty(mat_file)
    load(mat_file, 'OD', 'stim_mat');
    stim = StimClass('left');
    stim.SetData(stim_mat);
    stim.SetStates(stim_mat(:,[1 3]));

    dod.SetDataTimeSeries(OD);
end


dod = hmrR_MotionCorrectSplineSG(dod,mlActAuto,0.99,10,1);
dod = hmrR_BandpassFilt(dod,0,0.5);
[tIncAuto,tIncAutoCh] = hmrR_MotionArtifactByChannel(dod,probe,mlActMan,mlActAuto,tIncMan,0.5,1,20,5);
[stim,tRange] = hmrR_StimRejection(dod,stim,tIncAuto,tIncMan,cfg.trange);
dc = hmrR_OD2Conc(dod,probe,[6  6]);

[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,dA,tbasis,tHRF] = hmrR_GLM_time_basis(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,cfg.trange,cfg.AR_model,cfg.basis,cfg.gamma_parameters,cfg.rhoSD_ssThresh.homer,cfg.SS_mode,3);

Proc_data.dod = dod;
Proc_data.dA = dA;
Proc_data.time = dc.time;
Proc_data.tHRF = tHRF;
Proc_data.stim = stim;
Proc_data.tbasis = tbasis;
Proc_data.mlActAuto = mlActAuto;
Proc_data.dc = dc;
Proc_data.dcAvg = dcAvg;
Proc_data.probe = probe;