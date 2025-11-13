classdef SSDOTPipelineWindow < handle
    % Lightweight UI window that hosts the SSDOT pipeline actions.
    
    properties (Access = private)
        App     % matlab.ui.Figure or empty for standalone
        Fig     matlab.ui.Figure
        Grid    matlab.ui.container.GridLayout
        
        % Arrays for UI components
        Lamps   matlab.ui.control.Lamp
        
        % Buttons
        ButtonLoadAV        matlab.ui.control.Button
        ButtonDisplaySD     matlab.ui.control.Button
        ButtonLoadSens      matlab.ui.control.Button
        ButtonDisplaySensHead   matlab.ui.control.Button
        ButtonDisplaySensBrain   matlab.ui.control.Button
        ButtonPrepareOD     matlab.ui.control.Button
        ButtonGenSpatial    matlab.ui.control.Button
        ButtonDisplayKernels   matlab.ui.control.Button
        ButtonViewSingle    matlab.ui.control.Button
        ButtonRunRecon      matlab.ui.control.Button
        
        StatusLabel    matlab.ui.control.Label

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
            % obj.App = app;
            if ~isempty(app)
                obj.App = app;
            end
            obj.ssdot = ssdot;
            
            % Create child uifigure for the pipeline window
            obj.Fig = uifigure('Name','SSDOT Pipeline','Position',[200 150 700 420]);
            
            % Tie lifetime to main app (only if app exists)
            % if ~isempty(app) && isvalid(app)
            %     addlistener(app,'ObjectBeingDestroyed',@(~,~) obj.safeClose());
            % end
            
            % Layout: 6 rows x 3 columns (lamp, button1, button2)
            obj.Grid = uigridlayout(obj.Fig,[6 4]);
            obj.Grid.RowHeight    = {30,30,30,30,30,'1x'};
            obj.Grid.ColumnWidth  = {30,'1x','1x', '1x'};
            obj.Grid.Padding      = [12 12 12 12];
            obj.Grid.RowSpacing   = 8;
            obj.Grid.ColumnSpacing= 8;

            %% CREATE COMPONENT LIKE BUTTON, LISTBOX, ETC...
            % Create lamps for all 5 rows
            obj.Lamps = matlab.ui.control.Lamp.empty(0);
            for i = 1:5
                lamp = uilamp(obj.Grid);
                lamp.Layout.Row = i;
                lamp.Layout.Column = 1;
                lamp.Color = 'red';
                obj.Lamps(end+1) = lamp;
            end
            
            % --- Row 1: Buttons
            obj.ButtonLoadAV = uibutton(obj.Grid);
            obj.ButtonLoadAV.Text = 'Load AtlasViewer';
            obj.ButtonLoadAV.ButtonPushedFcn = @(s,e)obj.onLoadAV();
            obj.ButtonLoadAV.Layout.Row = 1;
            obj.ButtonLoadAV.Layout.Column = 2;
            
            obj.ButtonDisplaySD = uibutton(obj.Grid);
            obj.ButtonDisplaySD.Text = 'Display SD';
            obj.ButtonDisplaySD.ButtonPushedFcn = @(s,e)obj.onDisplaySD();
            obj.ButtonDisplaySD.Layout.Row = 1;
            obj.ButtonDisplaySD.Layout.Column = 3;
            
            % --- Row 2: Buttons
            obj.ButtonLoadSens = uibutton(obj.Grid);
            obj.ButtonLoadSens.Text = 'Load Sensitivity';
            obj.ButtonLoadSens.ButtonPushedFcn = @(s,e)obj.onLoadSensitivity();
            obj.ButtonLoadSens.Layout.Row = 2;
            obj.ButtonLoadSens.Layout.Column = 2;
            
            obj.ButtonDisplaySensHead = uibutton(obj.Grid);
            obj.ButtonDisplaySensHead.Text = 'Display Sensitivity Head';
            obj.ButtonDisplaySensHead.ButtonPushedFcn = @(s,e)obj.onDisplaySensitivityHead();
            obj.ButtonDisplaySensHead.Layout.Row = 2;
            obj.ButtonDisplaySensHead.Layout.Column = 3;

            obj.ButtonDisplaySensBrain = uibutton(obj.Grid);
            obj.ButtonDisplaySensBrain.Text = 'Display Sensitivity Brain';
            obj.ButtonDisplaySensBrain.ButtonPushedFcn = @(s,e)obj.onDisplaySensitivityBrain();
            obj.ButtonDisplaySensBrain.Layout.Row = 2;
            obj.ButtonDisplaySensBrain.Layout.Column = 4;
            
            % --- Row 3: Buttons       
            obj.ButtonGenSpatial = uibutton(obj.Grid);
            obj.ButtonGenSpatial.Text = 'Generate Spatial Base';
            obj.ButtonGenSpatial.ButtonPushedFcn = @(s,e)obj.onGenerateSpatialBase();
            obj.ButtonGenSpatial.Layout.Row = 3;
            obj.ButtonGenSpatial.Layout.Column = 2;

            obj.ButtonDisplayKernels = uibutton(obj.Grid);
            obj.ButtonDisplayKernels.Text = 'Display Spatial Base Centers';
            obj.ButtonDisplayKernels.ButtonPushedFcn = @(s,e)obj.onDisplaySpatialBaseCenters();
            obj.ButtonDisplayKernels.Layout.Row = 3;
            obj.ButtonDisplayKernels.Layout.Column = 3;
            
            % --- Row 4: Buttons      
            obj.ButtonViewSingle = uibutton(obj.Grid);
            obj.ButtonViewSingle.Text = 'View Single Time Base';
            obj.ButtonViewSingle.ButtonPushedFcn = @(s,e)obj.onViewSingleTimeBase();
            obj.ButtonViewSingle.Layout.Row = 4;
            obj.ButtonViewSingle.Layout.Column = 3;
            
            obj.ButtonPrepareOD = uibutton(obj.Grid);
            obj.ButtonPrepareOD.Text = 'Prepare OD Data';
            obj.ButtonPrepareOD.ButtonPushedFcn = @(s,e)obj.onPrepareOD();
            obj.ButtonPrepareOD.Layout.Row = 4;
            obj.ButtonPrepareOD.Layout.Column = 2;
            
            % --- Row 5:
            obj.ButtonRunRecon = uibutton(obj.Grid);
            obj.ButtonRunRecon.Text = 'Run Image Reconstruction';
            obj.ButtonRunRecon.ButtonPushedFcn = @(s,e)obj.onRunImageReconstruction();
            obj.ButtonRunRecon.Layout.Row = 5;
            obj.ButtonRunRecon.Layout.Column = [2 3];
            
            % --- Row 6: 
            obj.StatusLabel = uilabel(obj.Grid);
            obj.StatusLabel.Text = 'Ready';
            obj.StatusLabel.HorizontalAlignment = 'left';
            obj.StatusLabel.Layout.Row = 6;
            obj.StatusLabel.Layout.Column = [1 4];
            
            % Initial UI state
            obj.updateUI();
            
            % Debug output
            fprintf('Window created successfully!\n');
        end
        
        function safeClose(obj)
            if ~isempty(obj.Fig) && isvalid(obj.Fig)
                delete(obj.Fig);
            end
        end
        
        function updateUI(obj)
            % Gate buttons based on pipeline state
            obj.ButtonLoadAV.Enable = 'on';

            % 
            obj.canEnable(obj.ButtonDisplaySD, 'AtlasState');
            obj.canEnable(obj.ButtonLoadSens, 'AtlasState', 'selectedRun');
            obj.canEnable(obj.ButtonDisplaySensHead, 'AtlasState', 'Sensitivity_Matrix');
            obj.canEnable(obj.ButtonDisplaySensBrain, 'AtlasState', 'Sensitivity_Matrix');

            obj.canEnable(obj.ButtonGenSpatial, 'AtlasState', 'Sensitivity_Matrix');
            obj.canEnable(obj.ButtonDisplayKernels, 'AtlasState', 'Sensitivity_Matrix')
            obj.canEnable(obj.ButtonViewSingle, 'selectedRun');
            obj.canEnable(obj.ButtonPrepareOD, 'Sensitivity_Matrix', 'selectedRun', 'AtlasState', 'T');
            obj.canEnable(obj.ButtonRunRecon, 'A', 'G', 'Y', 'T', 'OD_SS', 'OD_drift', 'AtlasState')
            
        end

        function canEnable(obj, button, varargin)
            % Disable by default
            button.Enable = 'off';

            % Make sure ssdot is valid
            if isempty(obj.ssdot) || ~isvalid(obj.ssdot)
                return;
            end

            % Check each required property for non-empty value
            for i = 1:numel(varargin)
                fieldName = varargin{i};

                if ~isprop(obj.ssdot, fieldName)
                    return; % property missing
                end

                try
                    val = obj.ssdot.getVar(fieldName);
                catch
                    return; % getVar threw error (e.g., empty or invalid)
                end

                if isempty(val)
                    return; % property exists but not yet populated
                end
            end

            % All checks passed → enable button
            button.Enable = 'on';
            obj.setLampColor(button.Layout.Row, 'green');
        end
                
        %% CALL BACK FUNCTION
        function onLoadAV(obj)
            obj.StatusLabel.Text = 'Clicked Load AV';

            % Clear previous warnings
            lastwarn('');

            try
                LoggerWindow.log('info', 'Loading AtlasView...');

                % Run model function
                loadAV([], obj.ssdot);

                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg);   % forward model warning to logger
                end

                LoggerWindow.log('success', 'AtlasView loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
            end
        end
        
        function onDisplaySD(obj)
            try
                LoggerWindow.log('Info', 'Ploting SD');
                % run model
                plotSD(obj.App, obj.ssdot);
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'sd ploted');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
            end
            
        end
        
        function onLoadSensitivity(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                loadSens(obj.App, obj.ssdot);
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end
        
        function onDisplaySensitivityHead(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                plotSensitivity(obj.App, obj.ssdot, 'Head');
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end

        function onDisplaySensitivityBrain(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                plotSensitivity(obj.App, obj.ssdot, 'Brain');
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end
        
        
        function onGenerateSpatialBase(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                generate_spatial_base(obj.App, obj.ssdot);
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end
        

        function onDisplaySpatialBaseCenters(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                plot_kernels(obj.App, obj.ssdot)

                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end
        

        function onViewSingleTimeBase(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                plot_single_time_base(obj.App, obj.ssdot);
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end

        function onPrepareOD(obj)
            try
                LoggerWindow.log('info', 'Preparing/Processing Optical Densitiy');
                % run model
                prepare_OD_data(obj.App, obj.ssdot);

                [warnMsg, ~] = lastwarn;

                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end
                
                LoggerWindow.log('success', 'OD is ready');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
        end
        
        function onRunImageReconstruction(obj)
            try
                LoggerWindow.log('Info', 'Loading Sensitivity');
                % run model
                RunImageRecon(obj.App, obj.ssdot);
                % Detect model warnings
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    LoggerWindow.log('warning', warnMsg); 
                end

                LoggerWindow.log('success', 'Sensitivity is loaded');
                obj.updateUI();

            catch ME
                LoggerWindow.log('error', ME.message);
                uialert(obj.Fig, ME.message, 'Error');
                obj.StatusLabel.Text = 'Failed';
                
                rethrow(ME);
            end
            
        end

        
        function setStatus(obj, msg)
            if isvalid(obj.StatusLabel)
                obj.StatusLabel.Text = msg;
                fprintf('Status: %s\n', msg); % Debug output
            end
        end
        

        function setLampColor(obj, index, color)
            % Helper to change lamp color by index
            if index >= 1 && index <= length(obj.Lamps) && isvalid(obj.Lamps(index))
                obj.Lamps(index).Color = color;
                fprintf('Lamp %d changed to %s\n', index, color); % Debug output
            end
        end
    end


        methods (Access = public)
            function bringToFront(obj)
                if ~isempty(obj.Fig) && isvalid(obj.Fig)
                    figure(obj.Fig);
                end
            end
        end

end