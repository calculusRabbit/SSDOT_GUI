classdef LoadDisplayWindow < handle
    properties (Access = private)
        fig
        row = 1
        rowHeight = 30
    end

    methods
        function obj = LoadDisplayWindow()
            % create popup figure
            obj.fig = uifigure('Name', 'Load/Display');
        end

        function addStep(obj, description, actions)
            % description: text label for this row
            % actions: struct of name -> function_handle

            y = obj.fig.InnerPosition(4) - obj.row * obj.rowHeight;

            % Label on the left
            uilabel(obj.fig, 'Text', description, 'Position', [10, y, 120, 24]);

            % Create one button for each field in the struct
            x = 100;
            names = fieldnames(actions);
            for i = 1:numel(names)
                cb = actions.(names{i});
                uibutton(obj.fig, ...
                    'Text', names{i}, ...
                    'Position', [x, y, 100, 25], ...
                    'ButtonPushedFcn', @(~,~) cb());
                x = x + 100; % spacing between buttons
            end

            obj.row = obj.row + 1;
        end
    end
end
