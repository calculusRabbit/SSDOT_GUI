classdef ssdotClass < handle

    properties (Access = public)
        AtlasState
        Sensitivity_Matrix
        M
        A
        T
        G

        % Parameters/users' inputs
        cfg = struct()

        % Derived vectors/matrices for modeling
        Y = [] % stacked OD vector
        OD_SS = [] % matrix
        OD_drift = [] % matrix
        H
        b
        HTH
        HTY

        % container for all runs 
        RunList RunClass = RunClass.empty

        % selected run
        selectedRun RunClass = RunClass.empty


        meta struct = struct('created', datetime, 'notes', "")
    end

    methods (Access = public)
        function obj = ssdotClass()
            % default
            obj.cfg = struct( ...
                'threshold_brain', 5, ...
                'threshold_scalp', 5, ...
                'sigma_brain', 5, ...
                'sigma_scalp', 5, ...
                'threshold_sensitivity', -2);
        end

        %% set/save value function
        function setVar(obj, name, value)
            switch string(name)
                case "Run"
                    if isa(value, 'RunClass')
                        idx = value.iRun;
                        obj.RunList(idx) = value;
                        fprintf('Run asigned\n');
                    else
                        error('Expected a RunClass object for property "Run".');
                    end

                case 'selectedrun'
                    assert(isa(value,'RunClass'), 'selectedRun must be a RunClass object.');
                    obj.selectedRun = value;
                    
                case "AtlasState"
                    obj.checkAtlasState(value);
                    obj.AtlasState = value;
                    
                case "Sensitivity_Matrix"
                    obj.checkSentivityMatrix(value);
                    obj.Sensitivity_Matrix = value;

                case "M"
                    obj.checkM(value);
                    obj.M = value;

                case "A"
                    obj.checkA(value);
                    obj.A = value;

                case "T"
                    obj.checkT(value);      
                    obj.T = value;

                case "G"
                    obj.checkG(value);
                    obj.G = value;

                case "cfg"
                    obj.cfg = value;

                case "Y"
                    obj.checkY(Value);
                    obj.Y = value;

                case "OD_SS"
                    obj.checkOD_SS(value);
                    obj.OD_SS = value;

                case "OD_drift"
                    obj.checkOD_drift(value);
                    obj.OD_drift = value;

                case "meta"
                    obj.meta = value;s
                otherwise, error("Unknown var: %s", name);
            end
        end

        %% Get value 
        function out = getVar(obj, name)
            name = string(name);
            if ~isprop(obj, name)
                error('unknown field: %s', name);
            end

            if isempty(obj.(name))
                warning('field "%s" is empty.', name);
                out = [];
                return;
            end

            out = obj.(name);
        end

        
    end

    %% helper
    methods (Access = private)
        function checkAtlasState(~, AtlasState)
            assert(isstruct(AtlasState), 'AtlasState must be a struct (loaded from .mat).');
            assert(isfield(AtlasState,'probe'),   'AtlasState.probe missing.');
            assert(isfield(AtlasState,'fwmodel'), 'AtlasState.fwmodel missing.');
        end


        function checkSentivityMatrix(~, S) 
            %required field names
            requiredFields = ["Adot", "Adot_scalp", "E", "channels", "shortSepChLst", "Adot_orig", "Adot_scalp_orig"];

            % make sure sensitivity matrix is a struct
            if ~isstruct(S)
                LoggerWindow.log('error', 'Sensitivity_Matrix must be a struct.');
                return;
            end

            % check for required field for existence
            if ~all(isfield(S, requiredFields))
                LoggerWindow.log('error', 'Sensitivity_Matrix missing required fields.');
                return;
            end

            % Check for empty fields
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                if isempty(S.(fieldName))
                    LoggerWindow.log('warn', sprintf('Sensitivity_Matrix field "%s" is empty.', fieldName));;
                end
            end
        end


        function checkM(~, M)
            % Must be struct
            if ~isstruct(M)
                LoggerWindow.log('error', 'M must be a struct.');
                error('M must be a struct.');
            end
            

            % Required fields
            requiredFields = ["mask_brain", "mask_scalp"];
            if ~all(isfield(M, requiredFields))
                LoggerWindow.log('error', 'M missing required fields.');
                error('M missing required fields.');
            end

            % Check types and emptiness
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                val = M.(fieldName);

                if ~isnumeric(val)
                    LoggerWindow.log('error', sprintf('M.%s must be numeric.', fieldName));
                    error('M.%s must be numeric.', fieldName);
                end

                if isempty(val)
                    LoggerWindow.log('warn', sprintf('M.%s is empty.', fieldName));
                end
            end
        end


        function checkA(~, A)
            % Must be struct
            if ~isstruct(A)
                LoggerWindow.log('error', 'matrix A must be a struct.');
                error('A must be a struct.');
            end

            % Required top-level fields
            requiredMain = ["Amatrix_brain", "Amatrix_scalp"];
            if ~all(isfield(A, requiredMain))
                LoggerWindow.log('error', 'A missing required fields.');
                error('A missing required fields.');
            end

            % Required subfields for each part
            requiredSub = ["w1_HbO", "w1_HbR", "w2_HbO", "w2_HbR"];

            % Check Amatrix_brain
            brain = A.Amatrix_brain;
            if ~isstruct(brain)
                LoggerWindow.log('error', 'A.Amatrix_brain must be a struct.');
                error('A.Amatrix_brain must be a struct.');
            end
            if ~all(isfield(brain, requiredSub))
                LoggerWindow.log('error', 'A.Amatrix_brain missing required subfields.');
                error('A.Amatrix_brain missing required subfields.');
            end
            for i = 1:numel(requiredSub)
                fieldName = requiredSub(i);
                val = brain.(fieldName);
                if ~isnumeric(val)
                    LoggerWindow.log('error', sprintf('A.Amatrix_brain.%s must be numeric.', fieldName));
                    error('A.Amatrix_brain.%s must be numeric.', fieldName);
                end
                if isempty(val)
                    LoggerWindow.log('warn', sprintf('A.Amatrix_brain.%s is empty.', fieldName));
                end
            end

            % Check Amatrix_scalp
            scalp = A.Amatrix_scalp;
            if ~isstruct(scalp)
                LoggerWindow.log('error', 'A.Amatrix_scalp must be a struct.');
                error('A.Amatrix_scalp must be a struct.');
            end
            if ~all(isfield(scalp, requiredSub))
                LoggerWindow.log('error', 'A.Amatrix_scalp missing required subfields.');
                error('A.Amatrix_scalp missing required subfields.');
            end
            for i = 1:numel(requiredSub)
                fieldName = requiredSub(i);
                val = scalp.(fieldName);
                if ~isnumeric(val)
                    LoggerWindow.log('error', sprintf('A.Amatrix_scalp.%s must be numeric.', fieldName));
                    error('A.Amatrix_scalp.%s must be numeric.', fieldName);
                end
                if isempty(val)
                    LoggerWindow.log('warn', sprintf('A.Amatrix_scalp.%s is empty.', fieldName));
                end
            end
        end


        function checkG(~, G)
            if ~isstruct(G)
                LoggerWindow.log('error', 'G must be a struct.');
                error('G must be a struct.');
            end

            requiredFields = ["G_brain", "G_scalp"];
            if ~all(isfield(G, requiredFields))
                LoggerWindow.log('error', 'G is missing required fields.');
                error('G is missing required fields.');
            end

            % Each must be a numeric 2d matrix 
            for i = 1:numel(requiredFields)
                fname = requiredFields(i);
                val = G.(fname);

                if ~(isnumeric(val) && ismatrix(val))
                    LoggerWindow.log('error', sprintf('G.%s must be a numeric 2-D matrix.', fname));
                    error('G.%s must be a numeric 2-D matrix.', fname);
                end

                if isempty(val)
                    LoggerWindow.log('warn', sprintf('G.%s is empty.', fname));
                end
            end
        end


        function checkT(~, T)
             % Must be a struct
            if ~isstruct(T)
                LoggerWindow.log('error', 'T must be a struct.');
                error('T must be a struct.');
            end

            % Required fields
            requiredFields = ["T_HbO_brain", "T_HbR_brain", "t_HbO_brain", "t_HbR_brain", ...
                            "T_HbO_scalp", "T_HbR_scalp", "t_HbO_scalp", "t_HbR_scalp"];
            if ~all(isfield(T, requiredFields))
                LoggerWindow.log('error', 'T is missing one or more required fields.');
                error('T is missing one or more required fields.');
            end

            % Check each field
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                val = T.(fieldName);

                if ~isnumeric(val)
                    LoggerWindow.log('error', sprintf('T.%s must be numeric.', fieldName));
                    error('T.%s must be numeric.', fieldName);
                end

                if isempty(val)
                    LoggerWindow.log('warn', sprintf('T.%s is empty.', fieldName));
                end
            end
        end


        function checkY(~, Y)
            % Must be numeric and a vector
            if ~isnumeric(Y)
                LoggerWindow.log('error', 'Y must be numeric.');
                error('Y must be numeric.');
            end

            if ~isvector(Y)
                LoggerWindow.log('error', 'Y must be a vector.');
                error('Y must be a vector.');
            end

            % Basic emptiness check
            if isempty(Y)
                LoggerWindow.log('warn', 'Y is empty.');
                return;
            end

            % Optional sanity: finite values and not all zeros
            if any(~isfinite(Y))
                LoggerWindow.log('warn', 'Y contains NaN or Inf values.');
            end

            if all(Y == 0)
                LoggerWindow.log('warn', 'Y appears to be all zeros.');
            end
        end


        function checkOD_SS(~, OD_SS)
            if ~(isnumeric(OD_SS) && ismatrix(OD_SS))
                LoggerWindow.log('error', 'OD_SS must be a numeric 2-D matrix.');
                error('OD_SS must be a numeric 2-D matrix.');
            end

            if isempty(OD_SS)
                LoggerWindow.log('warn', 'OD_SS is empty.');
            end
        end


        function checkOD_drift(~, OD_drift)
            if ~(isnumeric(OD_drift) && ismatrix(OD_drift))
                LoggerWindow.log('error', 'OD_drift must be a numeric 2-D matrix.');
                error('OD_drift must be a numeric 2-D matrix.');
            end

            if isempty(OD_drift)
                LoggerWindow.log('warn', 'OD_drift is empty.');
            end
        end

        
    end
end


