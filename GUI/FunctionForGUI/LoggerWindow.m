classdef LoggerWindow < handle
    % Simple static logger window
    % Usage:
    %   LoggerWindow.log('info', "Hello");
    %   LoggerWindow.log('error', "Something failed");

    properties (Access = private)
        Fig
        Tbl
        Data = {}
        Styles
    end

    methods (Static)
        function log(level, msg)
            % Access the shared instance
            persistent inst
            if isempty(inst) || ~isvalid(inst)
                inst = LoggerWindow();
            end
            inst.addMessage(level, msg);
        end
    end

    methods (Access = private)
        function obj = LoggerWindow()
            % --- create UI ---
            obj.Fig = uifigure('Name', 'Logger', 'Position', [100 100 700 320]);
            obj.Tbl = uitable(obj.Fig, ...
                'Data', obj.Data, ...
                'ColumnName', {'Level','Message'}, ...
                'ColumnEditable', [false false], ...
                'ColumnWidth', {90,'1x'});
            obj.Tbl.FontName = 'Consolas';

            % --- color styles ---
            obj.Styles = struct( ...
                'info',    uistyle('FontColor',[0.15 0.15 0.15]), ... % gray
                'warn',    uistyle('FontColor',[0.85 0.55 0.00]), ... % orange
                'error',   uistyle('FontColor',[0.85 0.10 0.10]), ... % red
                'success', uistyle('FontColor',[0.10 0.55 0.10]) ...  % green
            );
        end

        function addMessage(obj, level, msg)
            % --- normalize ---
            level = lower(string(level));
            if ~isfield(obj.Styles, level)
                level = "info";
            end

            % --- append row ---
            newRow = {upper(char(level)), char(msg)};
            obj.Data = [obj.Data; newRow]; %#ok<AGROW>
            obj.Tbl.Data = obj.Data;

            % --- color style ---
            rowIdx = size(obj.Data,1);
            addStyle(obj.Tbl, obj.Styles.(level), 'row', rowIdx);
        end
    end
end
