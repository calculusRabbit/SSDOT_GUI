classdef LoggerWindow < handle
    %   LoggerWindow.log('info', 'Operation started');
    %   LoggerWindow.log('success', 'Operation completed');
    %   LoggerWindow.log('warning', 'Low memory');
    %   LoggerWindow.log('error', 'Failed to load file');
    %   LoggerWindow.show();    % Show the window
    %   LoggerWindow.hide();    % Hide the window
    %   LoggerWindow.clear();   % Clear all logs
    
    properties (Access = private)
        Fig
        TextArea
        SaveButton
        ClearButton
    end
    
    methods (Access = private)
        function obj = LoggerWindow()
            obj.createUI();
        end
        
        function createUI(obj)
            % Create the logger window
            obj.Fig = uifigure('Name', 'Log Window', ...
                              'Position', [100 100 700 450], ...
                              'Visible', 'off', ...
                              'CloseRequestFcn', @(~,~) obj.hideWindow());
            
            % main grid
            grid = uigridlayout(obj.Fig, [2 1]);
            grid.RowHeight = {'1x', 45};
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 10;
            
            % Text area for logs
            obj.TextArea = uitextarea(grid);
            obj.TextArea.Layout.Row = 1;
            obj.TextArea.Layout.Column = 1;
            obj.TextArea.Editable = 'off';
            obj.TextArea.FontName = 'Courier New';
            obj.TextArea.FontSize = 10;
            obj.TextArea.Value = {'Log window ready...'};
            
            % Button panel
            buttonPanel = uipanel(grid, 'BorderType', 'none');
            buttonPanel.Layout.Row = 2;
            buttonPanel.Layout.Column = 1;
            
            buttonGrid = uigridlayout(buttonPanel, [1 3]);
            buttonGrid.ColumnWidth = {100, 100, '1x'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.ColumnSpacing = 10;
            
            % Save Logs button
            obj.SaveButton = uibutton(buttonGrid, 'Text', 'Save Logs');
            obj.SaveButton.ButtonPushedFcn = @(~,~) obj.saveLogs();
            obj.SaveButton.Layout.Row = 1;
            obj.SaveButton.Layout.Column = 1;
            obj.SaveButton.Tooltip = 'Save all logs to a text file';
            
            % Clear button
            obj.ClearButton = uibutton(buttonGrid, 'Text', 'Clear Logs');
            obj.ClearButton.ButtonPushedFcn = @(~,~) obj.clearLogs();
            obj.ClearButton.Layout.Row = 1;
            obj.ClearButton.Layout.Column = 2;
            obj.ClearButton.Tooltip = 'Clear all log entries';
        end
        
        function addLogEntry(obj, level, message)
            % Add a log entry with timestamp and level
            timestamp = datetime('now', 'Format', 'HH:mm:ss');
            
            switch lower(level)
                case 'info'
                    prefix = '[INFO]   ';
                case 'success'
                    prefix = '[SUCCESS]';
                case 'warning'
                    prefix = '[WARNING]';
                case 'error'
                    prefix = '[ERROR]  ';
                otherwise
                    prefix = '[LOG]    ';
            end
            
            logEntry = sprintf('%s %s %s', string(timestamp), prefix, message);
            
            % Add to text area
            currentLogs = obj.TextArea.Value;
            currentLogs{end+1} = logEntry;
            obj.TextArea.Value = currentLogs;
            
            % Auto-scroll to bottom
            scroll(obj.TextArea, 'bottom');
            
            % Auto-show only for errors (user can still open manually anytime)
            if strcmp(obj.Fig.Visible, 'off') && strcmpi(level, 'error')
                obj.showWindow();
            end
            
            drawnow limitrate;
        end
        
        function saveLogs(obj)
            % Get current time as a datetime object
            dt = datetime('now'); 
            
            dateTimeString = string(dt, 'yyyy-MM-dd_HHmm'); 
            defaultName = sprintf('logs_%s.txt', dateTimeString);
            
            % ask user where to save
            [file, path] = uiputfile('*.txt', 'Save Logs As', defaultName);
            
            % if the user cancelled
            if file == 0, return; end
            
            logs = obj.TextArea.Value;
            fullPath = fullfile(path, file);
            
            % Write the log content to the selected file
            writecell(logs, fullPath, 'FileType', 'text', 'QuoteStrings', 'none');
            
            obj.addLogEntry('success', sprintf('Logs saved to: %s', fullPath));
        end
        
        function clearLogs(obj)
            % ask before clearing
            selection = uiconfirm(obj.Fig, ...
                'Are you sure you want to clear all log messages? This cannot be undone.', ...
                'Clear Logs', ...
                'Options', {'Clear', 'Cancel'});
            
            if ~strcmp(selection, 'Clear')
                return; 
            end

            obj.TextArea.Value = {'Log window cleared...'};
            obj.addLogEntry('info', 'Logs cleared');
        end
        
        function showWindow(obj)
            obj.Fig.Visible = 'on';
            figure(obj.Fig);
        end
        
        function hideWindow(obj)
            obj.Fig.Visible = 'off';
        end
        
        function visible = isWindowVisible(obj)
            visible = strcmp(obj.Fig.Visible, 'on');
        end
    end
    
    % ==== SINGLETON SUPPORT ====
    methods (Static, Access = private)
        function instance = shared(newInstance)
            persistent inst
            if nargin > 0
                inst = newInstance;
            end
            instance = inst;
        end
    end
    
    methods (Static)
        function log(level, message)
            instance = LoggerWindow.shared();
            if isempty(instance) || ~isvalid(instance)
                instance = LoggerWindow();
                LoggerWindow.shared(instance);
            end
            instance.addLogEntry(level, message);
        end
        
        function show()
            instance = LoggerWindow.shared();
            if isempty(instance) || ~isvalid(instance)
                instance = LoggerWindow();
                LoggerWindow.shared(instance);
            end
            instance.showWindow();
        end
        
        function hide()
            instance = LoggerWindow.shared();
            if ~isempty(instance) && isvalid(instance)
                instance.hideWindow();
            end
        end
        
        function clear()
            instance = LoggerWindow.shared();
            if ~isempty(instance) && isvalid(instance)
                instance.clearLogs();
            end
        end
        
        function visible = isVisible()
            instance = LoggerWindow.shared();
            visible = ~isempty(instance) && isvalid(instance) && instance.isWindowVisible();
        end
        
        function deleteInstance()
            instance = LoggerWindow.shared();
            if ~isempty(instance) && isvalid(instance)
                if isvalid(instance.Fig)
                    delete(instance.Fig);
                end
            end
            LoggerWindow.shared([]);
        end
    end
end