function simulate_pertubation(cfg)

%% Prepare savepath folder
resting_data_folder = cfg.restingpath;
cfg.savepath_rs = fullfile(cfg.savepath,'simulated_rs',num2str(cfg.scale_factor));
cfg.savepath_zeros = fullfile(cfg.savepath,'simulated_zeros',num2str(cfg.scale_factor));
cfg.savepath_snirf = fullfile(cfg.savepath,'simulated_snirf');
cfg.savepath_resting = fullfile(cfg.savepath,'resting');
mkdir(cfg.savepath_rs)
mkdir(cfg.savepath_zeros)
mkdir(cfg.savepath_snirf)
mkdir(cfg.savepath_resting)
mkdir(fullfile(cfg.savepath_zeros,'Subject'))
%% aquire brain model vertices
unilateral_data_path = cfg.unilateral_data_path;
At_file = 'atlasViewer.mat';
atlasViewer = load([unilateral_data_path, At_file]);
brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

%% make perturbation
axes_order = [2,1,3];
index_perturbation = [];

for i = 1:size(cfg.center,1)
    perturbation_center = cfg.center(i,:);
    distance = sqrt((brain_vertices(:,axes_order(1)) - perturbation_center(1)).^2 + ...
        (brain_vertices(:,axes_order(2)) - perturbation_center(2)).^2 + ...
        (brain_vertices(:,axes_order(3)) - perturbation_center(3)).^2);
    index_perturbation = [index_perturbation; find(distance<15)];
end

visualize_perturbation(cfg,atlasViewer,brain_vertices,index_perturbation,cfg.scale_factor)
%% read the resting state data
file = '2020-12-25_006_probe_correctSD.snirf';
acquired = SnirfClass([unilateral_data_path, file]);
rhoSD_ssThresh = 20;
ml     = acquired.data.GetMeasListSrcDetPairs();
SrcPos = acquired.probe.GetSrcPos();
DetPos = acquired.probe.GetDetPos();
lst = 1:size(ml,1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SrcPos(ml(lst(iML),1),:) - DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SrcPos(ml(lst(iML),1),:) + DetPos(ml(lst(iML),2),:)) / 2;
end
lstSS = lst(find(rhoSD<=rhoSD_ssThresh));
%% calculate OD
rs_folder = dir(fullfile(resting_data_folder,'*'));
subject_folder = setdiff({rs_folder([rs_folder.isdir]).name},{'.','..'});

for ii = 1:numel(subject_folder)
    mkdir([cfg.savepath_rs,filesep,subject_folder{ii}])
    mkdir([cfg.savepath_snirf,filesep,subject_folder{ii}])
    snirffile = dir(fullfile(resting_data_folder,subject_folder{ii},'*.snirf'));
    snirffilelist = {snirffile(~[snirffile.isdir]).name};
    for jj = 1:numel(snirffilelist)
        snirf = fullfile(resting_data_folder,subject_folder{ii},snirffilelist{jj});
        snirfObj = SnirfClass(snirf);
        fprintf('%s:\n',snirf)
        t = snirfObj.GetTime();
        % cut to the same length
        snirfObj.data.SetTime(t(1:2578))
        x = snirfObj.data.GetDataTimeSeries();
        snirfObj.data.SetDataTimeSeries(x(1:2578,:));

        snirfObj = rearrange_channels(snirfObj); % here should take the right channels out
        OD_rearrange = hmrR_Intensity2OD(snirfObj.data);
        OD_rearrange = OD_rearrange.GetDataTimeSeries;
        t_new = snirfObj.GetTime();
        dt = t_new(2) - t_new(1);
        nT = length(t_new);
        % time basis function
        tbasis = build_tbasis(dt);
        [T,stim] = build_Tbasis(nT,dt,tbasis);
        save(fullfile(cfg.savepath_resting, [subject_folder{ii},'.mat']),'OD_rearrange','T','stim','t_new')

        if ~isfile(fullfile(cfg.savepath_snirf, subject_folder{ii}, snirffilelist{jj}))
            snirfObj.probe.detectorPos2D(17,:) = [67.8, 57, 0];
            snirfObj.probe.detectorPos3D(17,:) = [67.8, 57, 0];
            measurementList_wl_1 = snirfObj.data.measurementList(1,50).copy;
            measurementList_wl_1.detectorIndex = 17;
            measurementList_wl_2 = snirfObj.data.measurementList(1,100).copy;
            measurementList_wl_2.detectorIndex = 17;

            measurementList_tmp(1, 1:51) = [snirfObj.data.measurementList(1,1:50), measurementList_wl_1];
            measurementList_tmp(1, 52:102) = [snirfObj.data.measurementList(1,51:100), measurementList_wl_2];
            snirfObj.data.measurementList = measurementList_tmp;

            intensity = snirfObj.data.GetDataTimeSeries();
            [~, col_1] = find(intensity(:,1:50)==max(intensity(:,1:50),[],'all'));
            [~, col_2] = find(intensity(:,51:100)==max(intensity(:,51:100),[],'all'));
            new_intensity = [intensity(:,1:50),...
                intensity(:,col_1),...
                intensity(:,51:100),...
                intensity(:,col_2)];
            snirfObj.data.SetDataTimeSeries(new_intensity);
            snirfObj.Save(fullfile(cfg.savepath_snirf, subject_folder{ii}, snirffilelist{jj}));
        end

        OD_SS_avg_w1 = mean(OD_rearrange(:,lstSS),2) * 2 / 19;
        OD_SS_avg_w2 = mean(OD_rearrange(:, 50 + lstSS),2) * 2 / 19;
        for tt = 1:length(T)
            dA = T{tt};
            mlActAuto{1,1} = ones(100,1);
            OD_HRF = generate_HRF(brain_vertices, scalp_vertices, index_perturbation,unilateral_data_path,dA,mlActAuto,cfg.scale_factor);
            OD_HRF = OD_HRF';
            OD = OD_rearrange + OD_HRF;
            OD = [OD(:, 1:50), ...
                OD_SS_avg_w1,...
                OD(:, 51:100),...
                OD_SS_avg_w2];
            
            stim_ts = stim{tt};
            stim_time = t(stim_ts == 1);
            stim_mat = zeros(length(stim_time),3);
            stim_mat(:,1) = stim_time;
            stim_mat(:,2) = 10*ones(length(stim_time),1);
            stim_mat(:,3) = ones(length(stim_time),1);
            
            save(fullfile(cfg.savepath_rs,subject_folder{ii},[snirffilelist{jj}(1:end-6) '_', int2str(tt),'.mat']), 'OD', 'stim_mat')
            
            OD = [OD_HRF(:,1:50), ...
                zeros(size(OD_SS_avg_w1)),...
                OD_HRF(:,51:100),...
                zeros(size(OD_SS_avg_w2))];
            save(fullfile(cfg.savepath_zeros,'Subject', ['Subject_', int2str(tt),'.mat']), 'OD', 'stim_mat')
        end
    end
end
end

function tbasis = build_tbasis(dt)
nB = 1; % one gamma basis
trange = [-2, 18];
nPre = round(trange(1)/dt);
nPost = round(trange(2)/dt);
tHRF = (1*nPre*dt:dt:nPost*dt)';
ntHRF = length(tHRF);
nConc = 2; % HbO and HbR
tbasis = zeros(ntHRF,nB,nConc);
paramsBasis = [0.1 3 7 0.1 3 7];
for iConc = 1:nConc
    tau = paramsBasis((iConc-1)*3+1);
    sigma = paramsBasis((iConc-1)*3+2);
    T = paramsBasis((iConc-1)*3+3);

    tbasis(:,1,iConc) = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
    lstNeg = find(tHRF<0);
    tbasis(lstNeg,1,iConc) = 0;

    if tHRF(1)<tau
        tbasis(1:round((tau-tHRF(1))/dt),1,iConc) = 0;
    end

    if T>0
        for ii=1:nB
            foo = conv(tbasis(:,ii,iConc),ones(round(T/dt),1)) / round(T/dt);
            tbasis(:,ii,iConc) = foo(1:ntHRF,1);
        end
    end
end
end

function [Amatrix,stim] = build_Tbasis(nT,dt,tbasis)

load('results/stim_sequences.mat','index','task_mat','rest_mat');
% for each resting state run, we are simulating three types of stimulus sequences
n = 3;
nCond = 1;% we are using only one type of stimulus
nB = 1; % one gamma basis

for i = 1:n
    task_array = task_mat(index(i),:);
    rest_array = rest_mat(index(i),:);
    % build onset time
    onset = zeros(nT,nCond);
    time_point = 5/dt;
    for j = 1:length(task_array)
        if j == 1
            onset(round(time_point)) = 1;
        else
            time_point = time_point + 5/dt + rest_array(j)/dt;
            onset(round(time_point)) = 1;
        end
    end
    % calculate dA
    dA=zeros(nT,nB*nCond,2);
    for iConc = 1:2
        iC = 0;
        for iCond=1:nCond
            for b=1:nB
                iC = iC + 1;
                if size(tbasis,3)==1
                    clmn = conv(onset(:,iCond),tbasis(:,b));
                else
                    clmn = conv(onset(:,iCond),tbasis(:,b,iConc));
                end
                clmn = clmn(1:nT);
                dA(:,iC,iConc) = clmn;
                onset = onset(1:nT);
            end
        end
    end
    Amatrix{i} = dA;
    stim{i} = onset;
end
end


function OD_HRF = generate_HRF(brain_vertices,scalp_vertices, index_perturbation,acquired_path,dA,mlActAuto, scale_factor)
%% building b matrix
g = zeros(size(brain_vertices,1),1);
g(index_perturbation) = 1;
% scale_factor = 1e-2; % for HbO, max activation = 0.4*scale_factor
b_brain_HbO = kron(g,squeeze(dA(:,:,1))'*scale_factor); % the max of dA is 0.4, so the max activation is 0.4e-3 here
b_brain_HbR = kron(g,squeeze(dA(:,:,2))'*(-0.5*scale_factor));
b_scalp_HbO = zeros(size(scalp_vertices,1),size(dA,1));
b_scalp_HbR = zeros(size(scalp_vertices,1),size(dA,1));

%% building A matrix
rhoSD_ssThresh = 0;
Sensitivity_Matrix = Get_A_dot(acquired_path,rhoSD_ssThresh,mlActAuto,0);
E = Sensitivity_Matrix.E;
Adot = Sensitivity_Matrix.Adot_orig;
Adot_scalp = Sensitivity_Matrix.Adot_scalp_orig;

Adot = double(Adot);
Adot_scalp = double(Adot_scalp);
A = Make_A_matrix(Adot, Adot_scalp, E,  []);

Amatrix_brain = A.Amatrix_brain;
Amatrix_scalp = A.Amatrix_scalp;

%% calculate x matrix
A_brain = [Amatrix_brain.w1_HbO,Amatrix_brain.w1_HbR;...
    Amatrix_brain.w2_HbO,Amatrix_brain.w2_HbR];
A_scalp = [Amatrix_scalp.w1_HbO,Amatrix_scalp.w1_HbR;...
    Amatrix_scalp.w2_HbO,Amatrix_scalp.w2_HbR];
OD_HRF = [A_brain, A_scalp] * [b_brain_HbO; b_brain_HbR; b_scalp_HbO; b_scalp_HbR ];

end

function visualize_perturbation(cfg,atlasViewer,brain_vertices,index_perturbation,scale_factor)

h = figure('name','perturbation_HbO');
axes_order = [2,1,3];
intensity = zeros(size(brain_vertices,1),1);
intensity(index_perturbation) = 0.4*scale_factor;
faces = atlasViewer.fwmodel.mesh.faces;
trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','interp','edgealpha',0, 'visible','on');
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
caxis([-1 1]*0.5*scale_factor)
title('HbO')
set(gcf,'position',[ 10         10         399         235])
axis image
if(~exist('light_onoff') | (exist('light_onoff') & strcmp(light_onoff,'on')))
    l = camlight;
    set(l,'Position',[50 2000 100]);

    l2 = camlight;
    set(l2,'Position',[50 -100 -100]);

    camlight(0,0);
end
lighting phong;
myColorMap = jet(256);
myColorMap(127:129,:) = 0.8;
colormap(myColorMap);
colorbar
hold on

savefig([cfg.savepath_rs,filesep,'pertubation_HbO.fig'])
close(h)

h = figure('name','perturbation_HbR');
axes_order = [2,1,3];
intensity = zeros(size(brain_vertices,1),1);
intensity(index_perturbation) = -0.4*0.5*scale_factor;
faces = atlasViewer.fwmodel.mesh.faces;
trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','interp','edgealpha',0, 'visible','on');
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
caxis([-1 1]*0.5*scale_factor)
title('HbR')
set(gcf,'position',[ 10         10         399         235])
axis image
if(~exist('light_onoff') | (exist('light_onoff') & strcmp(light_onoff,'on')))
    l = camlight;
    set(l,'Position',[50 2000 100]);

    l2 = camlight;
    set(l2,'Position',[50 -100 -100]);

    camlight(0,0);
end
lighting phong;
myColorMap = jet(256);
myColorMap(127:129,:) = 0.8;
colormap(myColorMap);
colorbar

savefig([cfg.savepath_rs,filesep,'pertubation_HbR.fig'])
close(h)

end