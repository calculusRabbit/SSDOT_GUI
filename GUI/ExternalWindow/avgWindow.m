classdef avgWindow < handle

    properties (Access = private)
        Fig
        MainGrid

        %Data
        Target
        Items

        %StatusPanel
        ProcessedLabel
        NotProcessedLabel
        ComputeAvgButton
        TargetsTable

    end


    methods (Static)
        function win = open(Target)
            win = avgWindow(Target);
        end
    end

    
    methods (Access = private)
        function obj = avgWindow(Target)
            if isempty(Target)
                error('AverageResultWindow:InvalidInput', 'Valid result requrie');
            end

            obj.Target = Target;
            obj.Items = Target.children;

            % determine level and get list to display

            
            obj.createWindow();
        end


        function createWindow(obj)
            % create main window
            obj.Fig = uifigure;
            obj.Fig.Name = 'Result';
            obj.Fig.Position = [100 100 820 610];

            % create grid layout
            obj.MainGrid = uigridlayout(obj.Fig, [2 1])
            obj.MainGrid.RowHeight = {'0.3x', '0.7x'};
            obj.MainGrid.RowSpacing = 10;

            obj.createStatusPannel();
            obj.updateStatusPannel();
            obj.createResultPannel();
        end


        function createStatusPannel(obj)
            %  Main Status Pannel
            StatusPanel = uipanel(obj.MainGrid, 'Title', 'Status');
            StatusPanel.Layout.Row = 1;
            StatusPanel.Layout.Column = 1;

            % Create a grid to split left and right panels
            Statusgrid = uigridlayout(StatusPanel, [1 2]);
            Statusgrid.ColumnWidth = {'0.5x', '0.5x'};
            Statusgrid.ColumnSpacing = 10;

            % Left Pannel
            leftGrid = uigridlayout(Statusgrid, [3 1]);
            leftGrid.Layout.Row = 1;
            leftGrid.Layout.Column = 1;

            obj.ProcessedLabel = uilabel(leftGrid);
            obj.ProcessedLabel.Layout.Row = 1;
            obj.ProcessedLabel.Layout.Column = 1;

            obj.NotProcessedLabel = uilabel(leftGrid);
            obj.NotProcessedLabel.Layout.Row = 2;
            obj.NotProcessedLabel.Layout.Column = 1;
            
            obj.ComputeAvgButton = uibutton(leftGrid);
            obj.ComputeAvgButton.Text = 'Calculate Average';
            obj.ComputeAvgButton.Layout.Row = 3;
            obj.ComputeAvgButton.Layout.Column = 1;

            % Right Pannel
            obj.TargetsTable = uitable(Statusgrid);
            obj.TargetsTable.Layout.Row = 1;
            obj.TargetsTable.Layout.Column = 2;
            obj.TargetsTable.ColumnName = {'Target', 'Status', 'Select'};
            obj.TargetsTable.ColumnEditable = [false false true];
        end


        function createResultPannel(obj)
            % Create result panel
            ResultsPanel = uipanel(obj.MainGrid, 'Title', 'Results');
            ResultsPanel.Layout.Row = 2;
            ResultsPanel.Layout.Column = 1;

            numItems = numel(obj.Items);

            % Grid that scrolls horizontally
            ResultsGrid = uigridlayout(ResultsPanel, [1 numItems]);
            ResultsGrid.ColumnWidth = repmat({250}, 1, numItems);
            ResultsGrid.RowHeight = {'1x'};
            ResultsGrid.Scrollable = 'on';
            ResultsGrid.ColumnSpacing = 10;

            % Add one plot per item
            for i = 1:numItems

                plotPanel = uipanel(ResultsGrid);
                plotPanel.Layout.Row = 1;
                plotPanel.Layout.Column = i;

                ax = uiaxes(plotPanel);
                ax.Position = [10 10 230 180];
                title(ax, obj.Target.children(i).name);

                % DEMO
                x = 1:10;
                y = rand(1,10);
                plot(ax, x, y, '-o');
                grid(ax, 'on');
            end
        end
    end


    methods(Access = private)

        function updateStatusPannel(obj)
            % update left side
            numItems = numel(obj.Target.children);
            processedCount = 0;
            statusList = strings(numItems, 1);

            for i = 1:numItems
                if ~isempty(obj.Target.children(i).ssdot.b) && ~isempty(obj.Target.children(i).ssdot)
                    statusList(i) = "Processed";
                    processedCount = processedCount + 1;
                else
                    statusList(i) = "Not Processed";
                end
            end
            notProcessedCount = numItems - processedCount;

            str1 = sprintf("✓ %d Processed", processedCount);
            str2 = sprintf("✗ %d Not Processed", notProcessedCount);

            obj.ProcessedLabel.Text = str1;
            obj.NotProcessedLabel.Text = str2;

            % update right side
            % table
            names = obj.getChildrenName();
            selected = false(numItems, 1);
            
            obj.TargetsTable.Data = [cellstr(names') , cellstr(statusList) , num2cell(selected)];

        end



        function output = getChildrenName(obj)
            childSize = numel(obj.Target.children);
            output = strings(1, childSize);

            switch lower(obj.Target.type)
                case 'sess'
                    nameIdx = 'iRun';
                case 'subj'
                    nameIdx = 'iSess';
                case 'group'
                    nameIdx = 'iSubj';
                otherwise
                    nameIdx = 'name';
            end

            for i = 1:childSize
                child = obj.Target.children(i);
                output(i) = string(child.type) + " " + string(child.(nameIdx));
            end
        end

    end
end