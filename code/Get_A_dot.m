function Sensitivity_Matrix = Get_A_dot(acquired_path,rhoSD_ssThresh,mlActAuto,spatially_regu)

Sensitivity_Matrix.Adot = [];
Sensitivity_Matrix.Adot_scalp = [];
Sensitivity_Matrix.E = [];
Sensitivity_Matrix.channels = [];
Sensitivity_Matrix.shortSepChLst = [];
Sensitivity_Matrix.Adot_orig = [];
Sensitivity_Matrix.Adot_scalp_orig = [];

% load At file
A_file          =   fullfile(acquired_path, 'fw', 'Adot.mat');
A_scalp_file    =   fullfile(acquired_path, 'fw', 'Adot_scalp.mat');
Adot        = load(A_file);
Adot_scalp  = load(A_scalp_file);
Adot = Adot.Adot;
Adot_scalp = Adot_scalp.Adot_scalp;

At_file = 'atlasViewer.mat';
atlasViewer = load([acquired_path,At_file]);
warning off

probe = atlasViewer.probe;
SD = convertProbe2SD(probe);

ml = SD.MeasList;
if size(ml,1) ~= size( mlActAuto{1,1},1)
    mlActAuto{1,1}([51 102],:) = [];
end
activeChLst = find(ml(:,4)==1 & mlActAuto{1,1}==1);

lst = find(ml(:,4)==1 & SD.MeasListAct==1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
end
% rhoSD_ssThresh = 0.8;
longSepChLst = lst(find(rhoSD>=rhoSD_ssThresh));
shortSepChLst = lst(find(rhoSD<rhoSD_ssThresh));
lstLS_all = [longSepChLst; longSepChLst+size(ml,1)/2]; % both wavelengths

if isempty(lstLS_all)
    menu(sprintf('All channels meet short separation threshold.\nYou need some long separation channels for image recon.\nPlease lower the threshold and retry.'), 'Okay');
    return;
end

[channels,~]=intersect(activeChLst,longSepChLst);
Adot_orig = Adot;
Adot_scalp_orig = Adot_scalp;

Adot = Adot(channels,:,:);
Adot_scalp = Adot_scalp(channels,:,:);


E = GetExtinctions([760 850]);
E = E/10; %convert from /cm to /mm  E raws: wavelength, columns 1:HbO, 2:HbR

Sensitivity_Matrix.E = E;
Sensitivity_Matrix.channels = channels;
Sensitivity_Matrix.shortSepChLst = shortSepChLst;
Sensitivity_Matrix.Adot_orig = Adot_orig;
Sensitivity_Matrix.Adot_scalp_orig = Adot_scalp_orig;

if spatially_regu
    A_file = fullfile(acquired_path,'fw', 'Adot_spatially_regu.mat');
    A_scalp_file = fullfile(acquired_path,'fw','Adot_scalp_spatially_regu.mat');
    ll_file = fullfile(acquired_path,'fw','ll.mat');
    ll_scalp_file = fullfile(acquired_path,'fw','ll_scalp.mat');
    if ~exist(A_file,'file') || ~exist(A_scalp_file,'file') || ~exist(ll_file,'file') || ~exist(ll_scalp_file,'file') 
        Adot_mat = Adot_orig;
        Adot_scalp_mat = Adot_scalp_orig;
        for i = 1:size(Adot_mat,3)
            A = squeeze(Adot_mat(:,:,i));
            ll_0=sum(A.^2,1); % diag(L) = diag(A'A) (shortcut)
            lambda2 = 0.1;
            ll=sqrt(ll_0+lambda2*max(ll_0)); % Adjust with Lambda2 cut-off value
            A=bsxfun(@rdivide,A,ll);
            Adot_mat(:,:,i) = A;
        end
        save(fullfile(acquired_path,'fw', 'll.mat'),'ll')
        for i = 1:size(Adot_scalp_mat,3)
            A = squeeze(Adot_scalp_mat(:,:,i));
            ll_0=sum(A.^2,1); % diag(L) = diag(A'A) (shortcut)
            lambda2 = 0.1;
            ll_scalp=sqrt(ll_0+lambda2*max(ll_0)); % Adjust with Lambda2 cut-off value
            A=bsxfun(@rdivide,A,ll_scalp);
            Adot_scalp_mat(:,:,i) = A;
        end
        save(fullfile(acquired_path,'fw', 'll_scalp.mat'),'ll_scalp')
        save(fullfile(acquired_path,'fw', 'Adot_spatially_regu.mat'),'Adot_mat')
        save(fullfile(acquired_path,'fw', 'Adot_scalp_spatially_regu.mat'),'Adot_scalp_mat')
    else
        load(fullfile(acquired_path,'fw', 'll.mat'),'ll')
        load(fullfile(acquired_path,'fw', 'll_scalp.mat'),'ll_scalp')
        load(fullfile(acquired_path,'fw', 'Adot_spatially_regu.mat'),'Adot_mat')
        load(fullfile(acquired_path,'fw', 'Adot_scalp_spatially_regu.mat'),'Adot_scalp_mat')
    end
    Adot = Adot_mat(channels,:,:);
    Adot_scalp= Adot_scalp_mat(channels,:,:);
    Sensitivity_Matrix.ll = ll;
    Sensitivity_Matrix.ll_scalp = ll_scalp;
end

Sensitivity_Matrix.Adot = Adot;
Sensitivity_Matrix.Adot_scalp = Adot_scalp;