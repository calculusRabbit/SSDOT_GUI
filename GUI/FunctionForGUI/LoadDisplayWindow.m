classdef LoadDisplayWindow < handle
    properties (Access = private)
        fig
        rowth = 1
    end

    methods
        function obj = LoadDisplayWindow()
            obj.fig = uifigure("Name","Load/Display");
        end

        function addStep(obj, description, loadFunc, displayFunc)
            maxHeight = obj.fig.InnerPosition(4);
            y = maxHeight - 30*obj.rowth; 

            % description
            uilabel(obj.fig, "Text", description, "Position", [10, y, 80, 22]);

            % Load button
            uibutton(obj.fig, "Text", "Load", ...
                "Position", [100, y, 60, 22], ...
                "ButtonPushedFcn", @(~,~) loadFunc());

            % Display button
            uibutton(obj.fig, "Text", "Display", ...
                "Position", [170, y, 70, 22], ...
                "ButtonPushedFcn", @(~,~) displayFunc());

            obj.rowth = obj.rowth + 1;
        end
    end
end
