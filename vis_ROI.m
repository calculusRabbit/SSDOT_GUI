%% original ROI
axes_order = [2,1,3];
file = fullfile('simulated_data', 'ROI', 'ROI', 'ROI.mat');
load(file)
index_nonROI_original = index_nonROI;
index_ROI_original = index_ROI;

av_data_path = ['data',filesep];
At_file = 'atlasViewer.mat';
atlasViewer = load([av_data_path, At_file]);
warning off
faces = atlasViewer.fwmodel.mesh.faces;
brain_vertices = atlasViewer.fwmodel.mesh.vertices; 

figure('name','brain');
intensity = zeros(size(brain_vertices,1), 1);
h = trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','black','facealpha',0.1,'edgealpha',0, 'visible','on');
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);

hold on
plot3(brain_vertices(index_ROI,axes_order(1)),brain_vertices(index_ROI,axes_order(2)),brain_vertices(index_ROI,axes_order(3)),'y.','linewidth',1,'markersize',13);
plot3(brain_vertices(index_nonROI,axes_order(1)),brain_vertices(index_nonROI,axes_order(2)),brain_vertices(index_nonROI,axes_order(3)),'c.','linewidth',1,'markersize',13);
% left
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])

set(gcf,'position',[10         10          313         242])
axis off
axis image

center1 = [94,177,135];
% plot3(center1(1),center1(2),center1(3),'rX','linewidth',1,'markersize',16);
%% ROI after process

file = fullfile('simulated_data', 'ROI_afterprocess', 'ROI', 'ROI.mat');
load(file)
index_nonROI_afterprocess = index_nonROI;
index_ROI_afterprocess = index_ROI;

figure('name','brain');
intensity = zeros(size(brain_vertices,1), 1);
h = trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','black','facealpha',0.1,'edgealpha',0, 'visible','on');
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
hold on
plot3(brain_vertices(index_ROI,axes_order(1)),brain_vertices(index_ROI,axes_order(2)),brain_vertices(index_ROI,axes_order(3)),'y.','linewidth',1,'markersize',13);
plot3(brain_vertices(index_nonROI,axes_order(1)),brain_vertices(index_nonROI,axes_order(2)),brain_vertices(index_nonROI,axes_order(3)),'c.','linewidth',1,'markersize',13);
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
set(gcf,'position',[10         10          313         242])
axis off
axis image

center_index = 17602;
center2 = brain_vertices(center_index,:);
% plot3(center2(axes_order(1)),center2(axes_order(2)),center2(axes_order(3)),'rX','linewidth',1,'markersize',16);

distance = sqrt((center1(1) - center2(axes_order(1))).^2 + ...
    (center1(2) - center2(axes_order(2))).^2 + ...
    (center1(3) - center2(axes_order(3))).^2)