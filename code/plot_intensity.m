function plot_intensity(cfg,faces, brain_vertices, int_at_pos, rotation)

axes_order = [2,1,3];

h = trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
      int_at_pos,'facecolor','interp','edgealpha',0, 'visible','on'); 
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
caxis(cfg.caxis_value)

if strcmp(rotation, 'L')
    view(90,0)
    camtarget([128.0, 132.0, 130.0])
    campos([128.0, 2238.8, 130.0])
    camup([-1.0, 0.0, 0.0])
elseif strcmp(rotation, 'R')
    view(-90,0)
    camtarget([128.0, 132.0, 130.0])
    campos([128.0, -2291.8, 130.0])
    camup([-1.0, 0.0, 0.0])
else
    fprintf('check the view positon\n')
end

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
% colorbar
axis image
axis off
hold on

if strcmp(cfg.loc_line,'on')
    perturbation_center = cfg.center;
    plot3([perturbation_center(1),perturbation_center(1)],[220,220],[40,220],'y-.','linewidth',1.5)
    plot3([60 200],[220 220],[perturbation_center(3),perturbation_center(3)],'y-.','linewidth',1.5)
end