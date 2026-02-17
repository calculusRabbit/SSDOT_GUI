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
        ResultAxes

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
            obj.updateResultsAxes();
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
            obj.ComputeAvgButton.ButtonPushedFcn = @(~,~) obj.onComputeAverage();

            % Right Pannel
            obj.TargetsTable = uitable(Statusgrid);
            obj.TargetsTable.Layout.Row = 1;
            obj.TargetsTable.Layout.Column = 2;
            obj.TargetsTable.ColumnName = {'Target', 'Status', 'Select'};
            obj.TargetsTable.ColumnEditable = [false false true];
        end


        function createResultPannel(obj)
            % create result panel to show each individual result that have been processed
            ResultsPannel = uipanel(obj.MainGrid, 'Title', 'Result');
            ResultsPannel.Layout.Row = 2;
            ResultsPannel.Layout.Column = 1;

            numItems = numel(obj.Target.children);

            % grid that scrolls horizontal
            ResultsGrid = uigridlayout(ResultsPannel, [1 numItems]);
            ResultsGrid.ColumnWidth = repmat({600}, 1, numItems);
            ResultsGrid.Scrollable = 'on';
            ResultsGrid.ColumnSpacing = 10;

            % add plot
            % Initialize cell array to hold axes for each child
            obj.ResultAxes = cell(numItems, 1);

            for i = 1:numItems
                % panel for each result
                individualPanel = uipanel(ResultsGrid);
                individualPanel.Layout.Row = 1;
                individualPanel.Layout.Column = i;
                individualPanel.Title = obj.Target.children(i).name;

                %
                grid = uigridlayout(individualPanel, [4 2]);
                grid.RowHeight = {300, 300, 300, 300};
                grid.Scrollable = 'on';
                grid.RowSpacing = 10;

                obj.ResultAxes{i} = gobjects(4, 2);
                for r = 1:4
                    for c = 1:2
                        % Create the axes
                        obj.ResultAxes{i}(r, c) = uiaxes(grid);
                        obj.ResultAxes{i}(r, c).Layout.Row = r;
                        obj.ResultAxes{i}(r, c).Layout.Column = c;
                        
                        % Set the axis limits behavior
                        axis(obj.ResultAxes{i}(r, c), 'tight');
                    end
                end
            end
        end

        function updateResultsAxes(obj)
            faces = [];
            brain_vertices = [];
            for i = 1:numel(obj.Target.children)
                if ~isempty(obj.Target.children(i).ssdot.AtlasState)
                    faces = obj.Target.children(i).ssdot.AtlasState.pialsurf.mesh_reduced.faces;
                    brain_vertices = obj.Target.children(i).ssdot.AtlasState.pialsurf.mesh_reduced.vertices;
                    break;
                end
            end

            caxis_value = [-1 1]*1e-9;

            for i = 1:numel(obj.Target.children)
                if isempty(obj.Target.children(i).ssdot.AtlasState)
                    continue;
                end

                HbO = obj.Target.children(i).ssdot.result.HbO;
                HbR = obj.Target.children(i).ssdot.result.HbR;

                for i_cond = 1:size(HbO, 2)
                    if i_cond > 2, break; end  % only 2 rows for HbO, 2 for HbR

                    intensity_HbO = HbO(:, i_cond);
                    intensity_HbR = HbR(:, i_cond);

                    % HbO in rows 1-2
                    obj.plot_intensity(obj.ResultAxes{i}(i_cond, 1), caxis_value, faces, brain_vertices, intensity_HbO, 'L');
                    title(obj.ResultAxes{i}(i_cond, 1), sprintf('HbO Cond %d - L', i_cond));
                    obj.plot_intensity(obj.ResultAxes{i}(i_cond, 2), caxis_value, faces, brain_vertices, intensity_HbO, 'R');
                    title(obj.ResultAxes{i}(i_cond, 2), sprintf('HbO Cond %d - R', i_cond));

                    % HbR in rows 3-4
                    obj.plot_intensity(obj.ResultAxes{i}(i_cond + 2, 1), caxis_value, faces, brain_vertices, intensity_HbR, 'L');
                    title(obj.ResultAxes{i}(i_cond + 2, 1), sprintf('HbR Cond %d - L', i_cond));
                    obj.plot_intensity(obj.ResultAxes{i}(i_cond + 2, 2), caxis_value, faces, brain_vertices, intensity_HbR, 'R');
                    title(obj.ResultAxes{i}(i_cond + 2, 2), sprintf('HbR Cond %d - R', i_cond));
                end
            end
        end
    end


    methods(Access = private)
        function onComputeAverage(obj)
            % Read checkmarks from Select column
            tableData    = obj.TargetsTable.Data;
            selectedMask = cell2mat(tableData(:, 3));
            selectedIdx  = find(selectedMask);

            if isempty(selectedIdx)
                uialert(obj.Fig, 'Please check at least one run.', 'No Selection');
                return;
            end

            % Accumulate HbO and HbR
            firstChild = obj.Target.children(selectedIdx(1));
            HbO_sum    = firstChild.ssdot.result.HbO;
            HbR_sum    = firstChild.ssdot.result.HbR;

            for k = 2:numel(selectedIdx)
                i       = selectedIdx(k);
                HbO_sum = HbO_sum + obj.Target.children(i).ssdot.result.HbO;
                HbR_sum = HbR_sum + obj.Target.children(i).ssdot.result.HbR;
            end

            HbO_avg = HbO_sum / numel(selectedIdx);
            HbR_avg = HbR_sum / numel(selectedIdx);

            % Grab mesh and plot using existing plot_Hb
            faces          = firstChild.ssdot.AtlasState.pialsurf.mesh_reduced.faces;
            brain_vertices = firstChild.ssdot.AtlasState.pialsurf.mesh_reduced.vertices;

            plot_Hb(obj.Fig, faces, brain_vertices, HbO_avg, HbR_avg);
        end

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


    methods (Static, Access = private)
        function plot_intensity(ax, caxis_value, faces, brain_vertices, int_at_pos, rotation)
            cla(ax);
            hold(ax, 'on');
            axes_order = [2, 1, 3];

            h = trisurf(faces, ...
                brain_vertices(:, axes_order(1)), ...
                brain_vertices(:, axes_order(2)), ...
                brain_vertices(:, axes_order(3)), ...
                int_at_pos, ...
                'FaceColor', 'interp', ...
                'EdgeAlpha', 0, ...
                'FaceAlpha', 1, ...
                'Parent', ax);
            set(h, 'DiffuseStrength', 0.9, 'SpecularStrength', 0.12, 'AmbientStrength', 0.2);

            if ~isempty(caxis_value)
                clim(ax, caxis_value); % clim is the modern replacement for caxis
            end

            if strcmp(rotation, 'L')
                view(ax, 90, 0);
                camtarget(ax, [128.0, 132.0, 130.0]);
                campos(ax, [128.0, 2238.8, 130.0]);
                camup(ax, [-1.0, 0.0, 0.0]);
            elseif strcmp(rotation, 'R')
                view(ax, -90, 0);
                camtarget(ax, [128.0, 132.0, 130.0]);
                campos(ax, [128.0, -2291.8, 130.0]);
                camup(ax, [-1.0, 0.0, 0.0]);
            end

            l = camlight(ax);
            set(l, 'Position', [50 2000 100]);
            l2 = camlight(ax);
            set(l2, 'Position', [50 -100 -100]);

            lighting(ax, 'phong');

            myColorMap = jet(256);
            myColorMap(127:129, :) = 0.8;
            colormap(ax, myColorMap);
            cb = colorbar(ax);
            cb.Position(4) = 0.5;
            cb.Position(2) = cb.Position(2) + 0.15;
            axis(ax, 'image');
            axis(ax, 'off');
        end
    end
end