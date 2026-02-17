classdef SSDOTPipelineWindow < handle
    
    properties (Access = private)
        App % this is for parent main.mlapp gui
        Fig matlab.ui.Figure
        MainGrid matlab.ui.container.GridLayout
        
        % property of ssdot class panels
        PropertyPanel matlab.ui.container.Panel
        PropertyGrid matlab.ui.container.GridLayout
        PropertyLamps
        
        % function group panels
        AtlasPanel matlab.ui.container.Panel
        SensitivityPanel matlab.ui.container.Panel
        SpatialPanel matlab.ui.container.Panel
        ODPanel matlab.ui.container.Panel
        ReconPanel matlab.ui.container.Panel
        
        % button
        ButtonLoadAV matlab.ui.control.Button
        ButtonDisplaySD matlab.ui.control.Button
        ButtonLoadSens matlab.ui.control.Button
        ButtonDisplaySensHead matlab.ui.control.Button
        ButtonDisplaySensBrain matlab.ui.control.Button
        ButtonGenSpatial matlab.ui.control.Button
        ButtonDisplayKernels matlab.ui.control.Button
        ButtonPrepareOD matlab.ui.control.Button
        ButtonViewSingle matlab.ui.control.Button
        ButtonRunRecon matlab.ui.control.Button
        
        %% User input:
        % for load Sensitivity
        FieldThresholdSensitivity matlab.ui.control.NumericEditField
        CheckSpatiallyRegu matlab.ui.control.CheckBox
        
        % generate Spatial Base
        FieldThresholdBrain matlab.ui.control.NumericEditField
        FieldThresholdScalp matlab.ui.control.NumericEditField
        FieldSigmaBrain matlab.ui.control.NumericEditField
        FieldSigmaScalp matlab.ui.control.NumericEditField
        
        % run Reconstruction
        DropdownDevice matlab.ui.control.DropDown
        DropdownRegularization matlab.ui.control.DropDown
        FieldAlpha matlab.ui.control.NumericEditField
        FieldBeta  matlab.ui.control.NumericEditField
        %%
        
        StatusLabel    matlab.ui.control.Label
        StatusLamp     matlab.ui.control.Lamp
        
    end


    properties (Access = public)
        ssdot
    end

    
    methods (Static)
        function h = open(app, ssdot)
            % Open with parent app
            h = SSDOTPipelineWindow(app, ssdot);
        end
        
        function h = openStandalone(ssdot)
            % Open without parent app (for testing)
            h = SSDOTPipelineWindow([], ssdot);
        end
    end
    
    methods (Access = private)
        function obj = SSDOTPipelineWindow(app, ssdot)
            % check ssdot object
            if isempty(ssdot) || ~isvalid(ssdot)
                error('SSDOTPipelineWindow:InvalidInput', 'Valid ssdot object required');
            end
            
            obj.ssdot = ssdot;
            
            if ~isempty(app)
                obj.App = app;
                addlistener(app, 'ObjectBeingDestroyed', @(src,evt) delete(obj.Fig)); %so if parent app close this one also close
            end
            
            % create main window
            obj.Fig = uifigure;
            obj.Fig.Name = 'SSDOT Pipeline';
            obj.Fig.Position = [100 100 1000 700];
            obj.Fig.Scrollable = 'on';
            obj.Fig.AutoResizeChildren = 'on';
            
            % main grid layout - responsive heights
            obj.MainGrid = uigridlayout(obj.Fig, [7 1]);
            obj.MainGrid.RowHeight = {'2x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            obj.MainGrid.Padding = [10 10 10 10];
            obj.MainGrid.RowSpacing = 8;
            
            % Create all panels for each sections
            obj.createPropertyStatusPanel();
            obj.createAtlasPanel();
            obj.createSensitivityPanel();
            obj.createSpatialPanel();
            obj.createODPanel();
            obj.createReconPanel();
            obj.createStatusBar();
            
            % initial UI state
            obj.updateUI();
        end


        %% CREATE UI
        function createPropertyStatusPanel(obj)
            % Property status panel - grid-based responsive layout
            obj.PropertyPanel = uipanel(obj.MainGrid, ...
                'Title', 'Data Status', ...
                'FontSize', 11, ...
                'FontWeight', 'bold');
            obj.PropertyPanel.Layout.Row = 1;
            obj.PropertyPanel.Layout.Column = 1;
            
            % Use grid layout for responsive design
            obj.PropertyGrid = uigridlayout(obj.PropertyPanel, [15 3]);
            obj.PropertyGrid.RowHeight = repmat({'fit'}, 1, 15);
            obj.PropertyGrid.ColumnWidth = {'fit', 'fit', '1x'};
            obj.PropertyGrid.Padding = [10 10 10 10];
            obj.PropertyGrid.RowSpacing = 4;
            obj.PropertyGrid.ColumnSpacing = 10;
            obj.PropertyGrid.Scrollable = 'on';
            
            % Property names and descriptions
            target = obj.ssdot.getVar('selectedRun');
            targetName = target.name;
            props = {
                'selectedRun', targetName;
                'AtlasState', 'AtlasViewer.mat';
                'Sensitivity_Matrix', 'Sensitivity data';
                'M', 'Masks';
                'A', 'A';
                'G', 'Spatial basis';
                'T', 'Temporal basis';
                'H', 'H matrix';
                'Y', 'Optical Density data vector';
                'OD_SS', 'Short separation Optical Density';
                'OD_drift', 'Drift Optical Density';
                'HTH', 'H^T * H';
                'HTY', 'H^T * Y';
                'Conc', 'Concentration';
                'b', 'Coefficients';
            };
            
            obj.PropertyLamps = struct();
            
            for i = 1:size(props, 1)
                % Lamp
                lamp = uilamp(obj.PropertyGrid);
                lamp.Color = [0.7 0.7 0.7];
                lamp.Layout.Row = i;
                lamp.Layout.Column = 1;
                obj.PropertyLamps.(props{i,1}) = lamp;

                % Property name label
                if i == 1
                    nameLabel = uilabel(obj.PropertyGrid, ...
                        'Text', target.type, ...
                        'FontWeight', 'bold');
                else
                    nameLabel = uilabel(obj.PropertyGrid, ...
                        'Text', props{i,1}, ...
                        'FontWeight', 'bold');
                end
                nameLabel.Layout.Row = i;
                nameLabel.Layout.Column = 2;
                
                % Description label
                descLabel = uilabel(obj.PropertyGrid, ...
                    'Text', props{i,2}, ...
                    'FontColor', [0.5 0.5 0.5]);
                descLabel.Layout.Row = i;
                descLabel.Layout.Column = 3;
            end   
        end
        
        function createAtlasPanel(obj)
            obj.AtlasPanel = uipanel(obj.MainGrid, 'Title', '1. AtlasViewer');
            obj.AtlasPanel.Layout.Row = 2;
            obj.AtlasPanel.Layout.Column = 1;
            
            grid = uigridlayout(obj.AtlasPanel, [1 3]);
            grid.ColumnWidth = {'fit', 'fit', '1x'};
            grid.Padding = [10 8 10 8];
            grid.ColumnSpacing = 8;
            
            obj.ButtonLoadAV = uibutton(grid, 'Text', 'Load Atlas');
            obj.ButtonLoadAV.ButtonPushedFcn = @(~,~) obj.onLoadAV();
            obj.ButtonLoadAV.Layout.Column = 1;
            
            obj.ButtonDisplaySD = uibutton(grid, 'Text', 'Display SD');
            obj.ButtonDisplaySD.ButtonPushedFcn = @(~,~) obj.onDisplaySD();
            obj.ButtonDisplaySD.Layout.Column = 2;
        end
        
        function createSensitivityPanel(obj)
            obj.SensitivityPanel = uipanel(obj.MainGrid, 'Title', '2. Sensitivity');
            obj.SensitivityPanel.Layout.Row = 3;
            obj.SensitivityPanel.Layout.Column = 1;
            
            grid = uigridlayout(obj.SensitivityPanel, [1 8]);
            grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 80, 'fit', 50, '1x'};
            grid.Padding = [10 8 10 8];
            grid.ColumnSpacing = 8;
            
            % Buttons
            obj.ButtonLoadSens = uibutton(grid, 'Text', 'Load Sens');
            obj.ButtonLoadSens.ButtonPushedFcn = @(~,~) obj.onLoadSensitivity();
            obj.ButtonLoadSens.Layout.Column = 1;
            
            obj.ButtonDisplaySensHead = uibutton(grid, 'Text', 'Display Head');
            obj.ButtonDisplaySensHead.ButtonPushedFcn = @(~,~) obj.onDisplaySensitivityHead();
            obj.ButtonDisplaySensHead.Layout.Column = 2;
            
            obj.ButtonDisplaySensBrain = uibutton(grid, 'Text', 'Display Brain');
            obj.ButtonDisplaySensBrain.ButtonPushedFcn = @(~,~) obj.onDisplaySensitivityBrain();
            obj.ButtonDisplaySensBrain.Layout.Column = 3;
            
            % Inline controls
            lbl1 = uilabel(grid, ...
                'Text', 'threshold_sens:', ...
                'HorizontalAlignment', 'right');
            lbl1.Layout.Column = 4;
            
            obj.FieldThresholdSensitivity = uieditfield(grid, 'numeric');
            obj.FieldThresholdSensitivity.Value = obj.ssdot.cfg.threshold_sensitivity;
            obj.FieldThresholdSensitivity.Limits = [-10 0];
            obj.FieldThresholdSensitivity.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldThresholdSensitivity.Layout.Column = 5;
            
            lbl2 = uilabel(grid, ...
                'Text', 'spatially_regu:', ...
                'HorizontalAlignment', 'right');
            lbl2.Layout.Column = 6;
            
            obj.CheckSpatiallyRegu = uicheckbox(grid, 'Text', '');
            obj.CheckSpatiallyRegu.Value = logical(obj.ssdot.cfg.spatially_regu);
            obj.CheckSpatiallyRegu.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.CheckSpatiallyRegu.Layout.Column = 7;
        end
        
        function createSpatialPanel(obj)
            obj.SpatialPanel = uipanel(obj.MainGrid, 'Title', '3. Spatial Base');
            obj.SpatialPanel.Layout.Row = 4;
            obj.SpatialPanel.Layout.Column = 1;
            
            grid = uigridlayout(obj.SpatialPanel, [1 10]);
            grid.ColumnWidth = {'fit', 'fit', 'fit', 75, 'fit', 75, 'fit', 75, 'fit', 75};
            grid.Padding = [10 8 10 8];
            grid.ColumnSpacing = 8;
            
            % Buttons
            obj.ButtonGenSpatial = uibutton(grid, 'Text', 'Generate');
            obj.ButtonGenSpatial.ButtonPushedFcn = @(~,~) obj.onGenerateSpatialBase();
            obj.ButtonGenSpatial.Layout.Column = 1;
            
            obj.ButtonDisplayKernels = uibutton(grid, 'Text', 'Show Centers');
            obj.ButtonDisplayKernels.ButtonPushedFcn = @(~,~) obj.onDisplaySpatialBaseCenters();
            obj.ButtonDisplayKernels.Layout.Column = 2;
            
            % User input
            lbl1 = uilabel(grid, ...
                'Text', 'threshold_brain:', ...
                'HorizontalAlignment', 'right');
            lbl1.Layout.Column = 3;
            
            obj.FieldThresholdBrain = uieditfield(grid, 'numeric');
            obj.FieldThresholdBrain.Value = obj.ssdot.cfg.threshold_brain;
            obj.FieldThresholdBrain.Limits = [0 Inf];
            obj.FieldThresholdBrain.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldThresholdBrain.Layout.Column = 4;
            
            lbl2 = uilabel(grid, ...
                'Text', 'threshold_scalp:', ...
                'HorizontalAlignment', 'right');
            lbl2.Layout.Column = 5;
            
            obj.FieldThresholdScalp = uieditfield(grid, 'numeric');
            obj.FieldThresholdScalp.Value = obj.ssdot.cfg.threshold_scalp;
            obj.FieldThresholdScalp.Limits = [0 Inf];
            obj.FieldThresholdScalp.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldThresholdScalp.Layout.Column = 6;
            
            lbl3 = uilabel(grid, ...
                'Text', 'sigma_brain:', ...
                'HorizontalAlignment', 'right');
            lbl3.Layout.Column = 7;
            
            obj.FieldSigmaBrain = uieditfield(grid, 'numeric');
            obj.FieldSigmaBrain.Value = obj.ssdot.cfg.sigma_brain;
            obj.FieldSigmaBrain.Limits = [0 Inf];
            obj.FieldSigmaBrain.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldSigmaBrain.Layout.Column = 8;
            
            lbl4 = uilabel(grid, ...
                'Text', 'sigma_scalp:', ...
                'HorizontalAlignment', 'right');
            lbl4.Layout.Column = 9;
            
            obj.FieldSigmaScalp = uieditfield(grid, 'numeric');
            obj.FieldSigmaScalp.Value = obj.ssdot.cfg.sigma_scalp;
            obj.FieldSigmaScalp.Limits = [0 Inf];
            obj.FieldSigmaScalp.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldSigmaScalp.Layout.Column = 10;
        end
        
        function createODPanel(obj)
            obj.ODPanel = uipanel(obj.MainGrid, 'Title', '4. OD Data');
            obj.ODPanel.Layout.Row = 5;
            obj.ODPanel.Layout.Column = 1;
            
            grid = uigridlayout(obj.ODPanel, [1 3]);
            grid.ColumnWidth = {'fit', 'fit', '1x'};
            grid.Padding = [10 8 10 8];
            grid.ColumnSpacing = 8;
            
            obj.ButtonPrepareOD = uibutton(grid, 'Text', 'Prepare OD');
            obj.ButtonPrepareOD.ButtonPushedFcn = @(~,~) obj.onPrepareOD();
            obj.ButtonPrepareOD.Layout.Column = 1;
            
            obj.ButtonViewSingle = uibutton(grid, 'Text', 'View Time Base');
            obj.ButtonViewSingle.ButtonPushedFcn = @(~,~) obj.onViewSingleTimeBase();
            obj.ButtonViewSingle.Layout.Column = 2;
        end
        
        function createReconPanel(obj)
            obj.ReconPanel = uipanel(obj.MainGrid, 'Title', '5. Reconstruction');
            obj.ReconPanel.Layout.Row = 6;
            obj.ReconPanel.Layout.Column = 1;
            
            grid = uigridlayout(obj.ReconPanel, [1 9]);
            grid.ColumnWidth = {'fit', 'fit', 75, 'fit', 50, 'fit', 75, 'fit', 75};
            grid.Padding = [10 8 10 8];
            grid.ColumnSpacing = 8;
            
            obj.ButtonRunRecon = uibutton(grid, 'Text', 'Run Recon');
            obj.ButtonRunRecon.ButtonPushedFcn = @(~,~) obj.onRunImageReconstruction();
            obj.ButtonRunRecon.Layout.Column = 1;
            
            % User input
            lbl1 = uilabel(grid, ...
                'Text', 'Device:', ...
                'HorizontalAlignment', 'right');
            lbl1.Layout.Column = 2;
            
            obj.DropdownDevice = uidropdown(grid);
            obj.DropdownDevice.Items = {'gpu', 'cpu'};
            obj.DropdownDevice.Value = obj.ssdot.cfg.device;
            obj.DropdownDevice.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.DropdownDevice.Layout.Column = 3;
            
            lbl2 = uilabel(grid, ...
                'Text', 'Reg:', ...
                'HorizontalAlignment', 'right');
            lbl2.Layout.Column = 4;
            
            obj.DropdownRegularization = uidropdown(grid);
            obj.DropdownRegularization.Items = {'1', '2', '3'};
            obj.DropdownRegularization.Value = num2str(obj.ssdot.cfg.regularization);
            obj.DropdownRegularization.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.DropdownRegularization.Layout.Column = 5;
            
            lbl3 = uilabel(grid, ...
                'Text', 'Alpha:', ...
                'HorizontalAlignment', 'right');
            lbl3.Layout.Column = 6;
            
            obj.FieldAlpha = uieditfield(grid, 'numeric');
            obj.FieldAlpha.Value = obj.ssdot.cfg.alpha;
            obj.FieldAlpha.Limits = [0 Inf];
            obj.FieldAlpha.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldAlpha.Layout.Column = 7;
            
            lbl4 = uilabel(grid, ...
                'Text', 'Beta:', ...
                'HorizontalAlignment', 'right');
            lbl4.Layout.Column = 8;
            
            obj.FieldBeta = uieditfield(grid, 'numeric');
            obj.FieldBeta.Value = obj.ssdot.cfg.beta;
            obj.FieldBeta.Limits = [0 Inf];
            obj.FieldBeta.ValueChangedFcn = @(~,~) obj.onConfigChanged();
            obj.FieldBeta.Layout.Column = 9;
        end
        
        function createStatusBar(obj)
            statusPanel = uipanel(obj.MainGrid, 'BorderType', 'line');
            statusPanel.Layout.Row = 7;
            statusPanel.Layout.Column = 1;
            
            grid = uigridlayout(statusPanel, [1 2]);
            grid.ColumnWidth = {'fit', '1x'};
            grid.Padding = [10 5 10 5];
            grid.ColumnSpacing = 8;
            
            obj.StatusLamp = uilamp(grid);
            obj.StatusLamp.Color = 'green';
            obj.StatusLamp.Layout.Column = 1;
            
            obj.StatusLabel = uilabel(grid);
            obj.StatusLabel.Text = 'Ready';
            obj.StatusLabel.FontSize = 10;
            obj.StatusLabel.FontWeight = 'bold';
            obj.StatusLabel.Layout.Column = 2;
        end
        
        %% 
        function onConfigChanged(obj)
            if isempty(obj.ssdot) || ~isvalid(obj.ssdot)
                return;
            end
            
            cfg = obj.ssdot.cfg;
            
            % Update from Sensitivity section
            cfg.threshold_sensitivity = obj.FieldThresholdSensitivity.Value;
            cfg.spatially_regu = double(obj.CheckSpatiallyRegu.Value);
            
            % Update from Spatial Base section
            cfg.threshold_brain = obj.FieldThresholdBrain.Value;
            cfg.threshold_scalp = obj.FieldThresholdScalp.Value;
            cfg.sigma_brain = obj.FieldSigmaBrain.Value;
            cfg.sigma_scalp = obj.FieldSigmaScalp.Value;
            
            % Update from Reconstruction section
            cfg.device = obj.DropdownDevice.Value;
            cfg.regularization = str2double(obj.DropdownRegularization.Value);
            cfg.alpha = obj.FieldAlpha.Value;
            cfg.beta = obj.FieldBeta.Value;
            
            % Update ssdot
            obj.ssdot.setVar('cfg', cfg);

            % debug
            fprintf('Config updated: device=%s, reg=%d, alpha=%.4f, beta=%.4f, thresh_sens=%.2f, spatial_regu=%d, th_brain=%.2f, th_scalp=%.2f, sigma_brain=%.2f, sigma_scalp=%.2f\n', ...
                cfg.device, cfg.regularization, cfg.alpha, cfg.beta, ...
                cfg.threshold_sensitivity, cfg.spatially_regu, ...
                cfg.threshold_brain, cfg.threshold_scalp, cfg.sigma_brain, cfg.sigma_scalp);
        end
        
        
        function updatePropertyLamps(obj)
            props = fieldnames(obj.PropertyLamps);
            
            for i = 1:length(props)
                propName = props{i};
                lamp = obj.PropertyLamps.(propName);
                
                try
                    if isprop(obj.ssdot, propName)
                        val = obj.ssdot.(propName);
                        if ~isempty(val)
                            lamp.Color = 'green';
                        else
                            lamp.Color = [0.7 0.7 0.7]; %gray
                        end
                    else
                        lamp.Color = [0.7 0.7 0.7];
                    end
                catch
                    lamp.Color = [0.7 0.7 0.7];
                end
            end
        end
        
        function canEnable(obj, button, varargin)
            button.Enable = 'off';
            
            if isempty(obj.ssdot) || ~isvalid(obj.ssdot)
                return;
            end
            
            for i = 1:numel(varargin)
                fieldName = varargin{i};
                
                if ~isprop(obj.ssdot, fieldName)
                    return;
                end
                
                try
                    val = obj.ssdot.getVar(fieldName);
                catch
                    return;
                end
                
                if isempty(val)
                    return;
                end
            end
            
            button.Enable = 'on';
        end

        function updateBtnTips(obj)
            % Tips for disabled buttons
            tips = {
                obj.ButtonDisplaySD, {'AtlasState'};
                obj.ButtonLoadSens, {'AtlasState', 'selectedRun'}
            };

            for i = 1:numel(tips(:,1))
                button = tips{i,1};
                reqFields = tips{i,2};

                if strcmp(button.Enable, 'off') && ~isempty(reqFields)
                    missing = {};
                    for j = 1:numel(reqFields)
                        if isempty(obj.ssdot.(reqFields{j}))
                            missing{end+1} = reqFields{j};
                        end
                    end
                    button.Tooltip = sprintf('Missing: %s\nComplete previous steps first.', strjoin(missing, ', '));
                end
            end
        end
        
        %%  Helper
        function executeWithLogging(obj, startMsg, func, successMsg)
            lastwarn('');
            obj.setStatus(startMsg, 'yellow');
            
            try
                LoggerWindow.log('info', startMsg);
                func();
                
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg);
                end
                
                LoggerWindow.log('success', successMsg);
                obj.setStatus(successMsg, 'green');
                obj.updateUI();
                
            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.setStatus('Operation failed', 'red');
                rethrow(ME);
            end
        end
        
        %% Callback functions
        function onLoadAV(obj)
            obj.executeWithLogging('Loading AtlasViewer...', ...
                @() loadAV([], obj.ssdot), ...
                'AtlasViewer loaded successfully');
        end
        
        function onDisplaySD(obj)
            obj.executeWithLogging('Plotting SD...', ...
                @() plotSD(obj.App, obj.ssdot), ...
                'SD plotted successfully');
        end
        
        function onLoadSensitivity(obj)
            obj.executeWithLogging('Loading sensitivity...', ...
                @() loadSens(obj.App, obj.ssdot), ...
                'Sensitivity loaded successfully');
        end
        
        function onDisplaySensitivityHead(obj)
            obj.executeWithLogging('Displaying sensitivity (Head)...', ...
                @() plotSensitivity(obj.App, obj.ssdot, 'Head'), ...
                'Sensitivity (Head) displayed successfully');
        end
        
        function onDisplaySensitivityBrain(obj)
            obj.executeWithLogging('Displaying sensitivity (Brain)...', ...
                @() plotSensitivity(obj.App, obj.ssdot, 'Brain'), ...
                'Sensitivity (Brain) displayed successfully');
        end
        
        function onGenerateSpatialBase(obj)
            obj.executeWithLogging('Generating spatial base...', ...
                @() generate_spatial_base(obj.App, obj.ssdot), ...
                'Spatial base generated successfully');
        end
        
        function onDisplaySpatialBaseCenters(obj)
            obj.executeWithLogging('Displaying spatial base centers...', ...
                @() plot_kernels(obj.App, obj.ssdot), ...
                'Spatial base centers displayed successfully');
        end
        
        function onViewSingleTimeBase(obj)
            obj.executeWithLogging('Viewing single time base...', ...
                @() plot_single_time_base(obj.App, obj.ssdot), ...
                'Time base displayed successfully');
        end
        
        function onPrepareOD(obj)
            obj.executeWithLogging('Preparing optical density data...', ...
                @() prepare_OD_data(obj.App, obj.ssdot), ...
                'OD data prepared successfully');
        end
        
        function onRunImageReconstruction(obj)
            obj.executeWithLogging('Running image reconstruction...', ...
                @() RunImageRecon(obj.App, obj.ssdot), ...
                'Image reconstruction completed successfully');
        end
        
        function setStatus(obj, msg, color)
            if nargin < 3
                color = 'green';
            end
            
            if isvalid(obj.StatusLabel)
                obj.StatusLabel.Text = msg;
                obj.StatusLamp.Color = color;
                drawnow;
            end
        end
    end
    
    methods (Access = public)
        function bringToFront(obj)
            if ~isempty(obj.Fig) && isvalid(obj.Fig)
                figure(obj.Fig);
            end
        end

        function updateUI(obj)
            % Update button states
            obj.canEnable(obj.ButtonLoadAV);
            obj.canEnable(obj.ButtonDisplaySD, 'AtlasState');
            obj.canEnable(obj.ButtonLoadSens, 'AtlasState', 'selectedRun');
            obj.canEnable(obj.ButtonDisplaySensHead, 'AtlasState', 'Sensitivity_Matrix');
            obj.canEnable(obj.ButtonDisplaySensBrain, 'AtlasState', 'Sensitivity_Matrix');
            obj.canEnable(obj.ButtonGenSpatial, 'AtlasState', 'Sensitivity_Matrix');
            obj.canEnable(obj.ButtonDisplayKernels, 'AtlasState', 'Sensitivity_Matrix', 'G');
            obj.canEnable(obj.ButtonViewSingle, 'selectedRun');
            obj.canEnable(obj.ButtonPrepareOD, 'Sensitivity_Matrix', 'selectedRun', 'AtlasState', 'T');
            obj.canEnable(obj.ButtonRunRecon, 'A', 'G', 'Y', 'T', 'OD_SS', 'OD_drift', 'AtlasState');
            
            obj.updatePropertyLamps();
            obj.updateBtnTips();
        end

        function valid = isWindowValid(obj)
            % Check if window is still open
            valid = ~isempty(obj.Fig) && isvalid(obj.Fig);
        end

        function close(obj)
            % Close the pipeline window
            if ~isempty(obj.Fig) && isvalid(obj.Fig)
                delete(obj.Fig);
            end
        end
    end
end