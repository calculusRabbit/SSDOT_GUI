classdef LoggerWindow < handle

    %   LoggerWindow.log('info', 'Loading data...');
    %   LoggerWindow.log('success', 'Operation complete!');
    %   LoggerWindow.log('error', 'Something went wrong');
    %   LoggerWindow.log('warning', 'Low memory detected');
    %   LoggerWindow.log('debug', 'Debug information');
    
    properties (Access = private)
        Figure
        TextArea
        ClearButton
        SaveButton
    end
    
    properties (Access = private, Constant)
        MaxLines = 1000  % Maximum number of log lines to keep
    end
    
    methods (Access = private)
        function obj = LoggerWindow()
            % Private constructor - enforces singleton pattern
            obj.createUI();
        end
    end
    
    methods (Access = public)
        function delete(obj)
            % Destructor - cleanup
            if isvalid(obj.Figure)
                delete(obj.Figure);
            end
        end
        
        function logMessage(obj, level, message)
            % Add a log message to the window
            
            % Get timestamp
            timestamp = datetime("now");
            
            % Format level string with icon
            switch lower(level)
                case 'info'
                    levelStr = '[INFO]';
                    icon = '▸';
                case 'success'
                    levelStr = '[OK]';
                    icon = '✓';
                case 'warning'
                    levelStr = '[WARN]';
                    icon = '⚠';
                case 'error'
                    levelStr = '[ERROR]';
                    icon = '✖';
                case 'debug'
                    levelStr = '[DEBUG]';
                    icon = '◦';
                otherwise
                    levelStr = '[LOG]';
                    icon = '•';
            end
            
            % Create log entry
            logEntry = sprintf('[%s] %s %-7s %s', timestamp, icon, levelStr, message);
            
            % Print to MATLAB console
            fprintf('%s\n', logEntry);
            
            % Add to text area if window exists
            if isvalid(obj.Figure) && isvalid(obj.TextArea)
                currentLog = obj.TextArea.Value;
                
                % Initialize or append
                if isempty(currentLog) || (numel(currentLog) == 1 && isempty(currentLog{1}))
                    obj.TextArea.Value = {logEntry};
                else
                    % Limit number of lines (prevent memory issues)
                    if numel(currentLog) >= obj.MaxLines
                        currentLog = currentLog(end - obj.MaxLines + 2 : end);
                    end
                    obj.TextArea.Value = [currentLog; {logEntry}];
                end
                
                drawnow;
            end
        end
        
        function clearLog(obj)
            % Clear all log messages
            if isvalid(obj.TextArea)
                obj.TextArea.Value = {''};
                fprintf('Logger cleared\n');
            end
        end
    end
    
    methods (Access = private)
        function createUI(obj)
            % Create the logger UI window
            
            obj.Figure = uifigure('Name', 'SSDOT Logger', ...
                                  'Position', [700 100 650 550]);
            
            % Main grid layout
            grid = uigridlayout(obj.Figure, [3 2]);
            grid.RowHeight = {'fit', '1x', 'fit'};
            grid.ColumnWidth = {'1x', '1x'};
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 8;
            
            % Title
            titleLabel = uilabel(grid, ...
                                 'Text', 'SSDOT Pipeline Logger', ...
                                 'FontSize', 14, ...
                                 'FontWeight', 'bold', ...
                                 'HorizontalAlignment', 'center');
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = [1 2];
            
            % Text area for logs
            obj.TextArea = uitextarea(grid, ...
                                      'Value', {''}, ...
                                      'Editable', 'off', ...
                                      'FontName', 'Courier New', ...
                                      'FontSize', 10);
            obj.TextArea.Layout.Row = 2;
            obj.TextArea.Layout.Column = [1 2];
            
            % Clear button
            obj.ClearButton = uibutton(grid, ...
                                       'Text', 'Clear Log', ...
                                       'ButtonPushedFcn', @(~,~) obj.clearLog());
            obj.ClearButton.Layout.Row = 3;
            obj.ClearButton.Layout.Column = 1;
            
            % Save button
            obj.SaveButton = uibutton(grid, ...
                                      'Text', 'Save Log to File', ...
                                      'ButtonPushedFcn', @(~,~) obj.saveLog());
            obj.SaveButton.Layout.Row = 3;
            obj.SaveButton.Layout.Column = 2;
        end
        
        function saveLog(obj)
            % Save log contents to a text file
            try
                [file, path] = uiputfile('*.txt', 'Save Log File', ...
                                         ['ssdot_log_' datestr(now, 'yyyymmdd_HHMMSS') '.txt']);
                if file ~= 0
                    fullPath = fullfile(path, file);
                    logContent = obj.TextArea.Value;
                    
                    % Write to file
                    fid = fopen(fullPath, 'w');
                    if fid == -1
                        error('Could not open file for writing');
                    end
                    
                    for i = 1:numel(logContent)
                        fprintf(fid, '%s\n', logContent{i});
                    end
                    fclose(fid);
                    
                    obj.logMessage('success', ['Log saved to: ' fullPath]);
                end
            catch ME
                uialert(obj.Figure, ME.message, 'Save Error');
                obj.logMessage('error', ['Failed to save log: ' ME.message]);
            end
        end
    end
    
    %% Static Methods (Singleton Pattern)
    methods (Static)
        function obj = getInstance()
            % Get or create the singleton logger instance
            persistent instance
            
            if isempty(instance) || ~isvalid(instance)
                instance = LoggerWindow();
            end
            
            obj = instance;
        end
        
        function log(level, message)         
            logger = LoggerWindow.getInstance();
            logger.logMessage(level, message);
        end
        
        function clear()
            % Clear all log messages
            % Usage: LoggerWindow.clear();
            logger = LoggerWindow.getInstance();
            logger.clearLog();
        end
    end
end