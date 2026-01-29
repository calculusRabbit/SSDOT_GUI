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
        
        % these data from runImageRecon
        HTH
        HTY
        Conc
        b

        % container for all runs 
        RunList RunClass = RunClass.empty

        % selected run
        selectedRun 

        % level: run/sess/subj/group
        level

        meta struct = struct('created', datetime, 'notes', "")
    end

    methods (Access = public)
        function obj = ssdotClass()
            % Constructor - initialize default configuration
            addpath(genpath('code'));
            obj.cfg = struct( ...
                'threshold_brain', 5, ...
                'threshold_scalp', 5, ...
                'sigma_brain', 5, ...
                'sigma_scalp', 5, ...
                'spatially_regu', 0, ...
                'threshold_sensitivity', -2, ...
                'alpha', 1e-2,...
                'beta', 0, ...
                'regularization', 1, ...
                'device', 'gpu');
            
            fprintf('ssdotClass initialized\n');
        end

        %% Set/save value function
        function setVar(obj, name, value)
            switch string(name)
                case "Run"
                    obj.checkRun(value);
                    idx = value.iRun;
                    obj.RunList(idx) = value;
                    
                case 'selectedrun'
                    obj.checkSelectedRun(value);
                    obj.selectedRun = value;

                case 'selectedRun'
                    obj.checkSelectedRun(value);
                    obj.selectedRun = value;

                case 'level'
                    obj.checkLevel(value)
                    obj.level = string(value);
                    
                case "AtlasState"
                    obj.checkAtlasState(value);
                    obj.AtlasState = value;
                    fprintf('AtlasState loaded successfully\n');
                    
                case "Sensitivity_Matrix"
                    obj.checkSentivityMatrix(value);
                    obj.Sensitivity_Matrix = value;
                    fprintf('Sensitivity_Matrix loaded successfully\n');

                case "M"
                    obj.checkM(value);
                    obj.M = value;
                    fprintf('M (mask) set successfully\n');

                case "A"
                    obj.checkA(value);
                    obj.A = value;
                    fprintf('A (sensitivity matrix) set successfully\n');

                case "T"
                    obj.checkT(value);      
                    obj.T = value;
                    fprintf('T (temporal basis) set successfully\n');

                case "G"
                    obj.checkG(value);
                    obj.G = value;
                    fprintf('G (spatial basis) set successfully\n');
                
                case "H"
                    obj.H = value;
                    fprintf('H matrix is set\n');

                case "cfg"
                    obj.cfg = value;
                    fprintf('Configuration updated\n');

                case "Y"
                    obj.checkY(value);
                    obj.Y = value;
                    fprintf('Y (data vector) set successfully\n');

                case "OD_SS"
                    obj.checkOD_SS(value);
                    obj.OD_SS = value;
                    fprintf('OD_SS set successfully\n');

                case "OD_drift"
                    obj.checkOD_drift(value);
                    obj.OD_drift = value;
                    fprintf('OD_drift set successfully\n');

                case "HTH"
                    obj.HTH = value;
                
                case "HTY"
                    obj.HTY = value;

                case "Conc"
                    obj.Conc = value;

                case "b"
                    obj.b = value;

                case "meta"
                    obj.meta = value;
                    fprintf('Metadata updated\n');
                    
                otherwise
                    error('SSDOT:UnknownVariable', 'Unknown variable: %s', name);
            end
        end

        %% Get value 
        function out = getVar(obj, name)
            name = string(name);
            
            if ~isprop(obj, name)
                error('SSDOT:UnknownField', 'Unknown field: %s', name);
            end

            if isempty(obj.(name))
                warning('SSDOT:EmptyField', 'Field "%s" is empty', name);
                out = [];
                return;
            end

            out = obj.(name);
        end

        % Save instance to .mat file
        function saveToFile(obj, filename)
            % optional: update metadata
            obj.meta.created = datetime("now");

            ssdot = obj;  
            save(filename, 'ssdot');
        end
    end

    %% Load functions
    methods (Access = public)
    end

    %% Helper validation methods
    methods (Access = private)
        function checkLevel(~, level)
            
        end

        function checkAtlasState(~, AtlasState)
            if ~isstruct(AtlasState)
                error('SSDOT:InvalidAtlasState', ...
                      'AtlasState must be a struct (loaded from .mat)');
            end
            
            if ~isfield(AtlasState, 'probe')
                error('SSDOT:MissingField', 'AtlasState.probe missing');
            end
            
            if ~isfield(AtlasState, 'fwmodel')
                error('SSDOT:MissingField', 'AtlasState.fwmodel missing');
            end
        end

        function checkRun(~, run)
            if ~isa(run, 'RunClass')
                error('SSDOT:InvalidType', ...
                      'Expected a RunClass object for property "Run"');
            end
        end

        function checkSelectedRun(~, node) 
            validClasses = {'RunClass', 'GroupClass', 'SubjClass', 'SessClass'}
            isValid = false;
            for i = 1:numel(validClasses)
                if isa(node, validClasses{i})
                     isValid = true;
                end
            end

            if ~isValid
                error('SSDOT:InvalidType', ...
                    'Expected RunClass, GroupClass, SubjClass, or SessClass');
            end
        end

        function checkSentivityMatrix(~, S) 
            % Required field names
            requiredFields = ["Adot", "Adot_scalp", "E", "channels", ...
                             "shortSepChLst", "Adot_orig", "Adot_scalp_orig"];

            if ~isstruct(S)
                error('SSDOT:InvalidType', 'Sensitivity_Matrix must be a struct');
            end

            % Check for required fields
            missingFields = requiredFields(~isfield(S, requiredFields));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'Sensitivity_Matrix missing required fields: %s', ...
                      strjoin(missingFields, ', '));
            end

            % Warn about empty fields
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                if isempty(S.(fieldName))
                    warning('SSDOT:EmptyField', ...
                            'Sensitivity_Matrix field "%s" is empty', fieldName);
                end
            end
        end

        function checkM(~, M)
            if ~isstruct(M)
                error('SSDOT:InvalidType', 'M must be a struct');
            end

            % Required fields
            requiredFields = ["mask_brain", "mask_scalp"];
            missingFields = requiredFields(~isfield(M, requiredFields));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'M missing required fields: %s', ...
                      strjoin(missingFields, ', '));
            end

            % Check types and emptiness
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                val = M.(fieldName);

                if ~isnumeric(val)
                    error('SSDOT:InvalidType', 'M.%s must be numeric', fieldName);
                end

                if isempty(val)
                    warning('SSDOT:EmptyField', 'M.%s is empty', fieldName);
                end
            end
        end

        function checkA(~, A)
            if ~isstruct(A)
                error('SSDOT:InvalidType', 'Matrix A must be a struct');
            end

            % Required top-level fields
            requiredMain = ["Amatrix_brain", "Amatrix_scalp"];
            if ~all(isfield(A, requiredMain))
                error('SSDOT:MissingFields', 'A missing required fields');
            end

            % Required subfields for each part
            requiredSub = ["w1_HbO", "w1_HbR", "w2_HbO", "w2_HbR"];

            % Check Amatrix_brain
            brain = A.Amatrix_brain;
            if ~isstruct(brain)
                error('SSDOT:InvalidType', 'A.Amatrix_brain must be a struct');
            end
            
            missingFields = requiredSub(~isfield(brain, requiredSub));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'A.Amatrix_brain missing subfields: %s', ...
                      strjoin(missingFields, ', '));
            end
            
            for i = 1:numel(requiredSub)
                fieldName = requiredSub(i);
                val = brain.(fieldName);
                if ~isnumeric(val)
                    error('SSDOT:InvalidType', ...
                          'A.Amatrix_brain.%s must be numeric', fieldName);
                end
                if isempty(val)
                    warning('SSDOT:EmptyField', ...
                            'A.Amatrix_brain.%s is empty', fieldName);
                end
            end

            % Check Amatrix_scalp
            scalp = A.Amatrix_scalp;
            if ~isstruct(scalp)
                error('SSDOT:InvalidType', 'A.Amatrix_scalp must be a struct');
            end
            
            missingFields = requiredSub(~isfield(scalp, requiredSub));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'A.Amatrix_scalp missing subfields: %s', ...
                      strjoin(missingFields, ', '));
            end
            
            for i = 1:numel(requiredSub)
                fieldName = requiredSub(i);
                val = scalp.(fieldName);
                if ~isnumeric(val)
                    error('SSDOT:InvalidType', ...
                          'A.Amatrix_scalp.%s must be numeric', fieldName);
                end
                if isempty(val)
                    warning('SSDOT:EmptyField', ...
                            'A.Amatrix_scalp.%s is empty', fieldName);
                end
            end
        end

        function checkG(~, G)
            if ~isstruct(G)
                error('SSDOT:InvalidType', 'G must be a struct');
            end

            requiredFields = ["G_brain", "G_scalp"];
            missingFields = requiredFields(~isfield(G, requiredFields));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'G missing required fields: %s', ...
                      strjoin(missingFields, ', '));
            end

            % Each must be a numeric 2D matrix 
            for i = 1:numel(requiredFields)
                fname = requiredFields(i);
                val = G.(fname);

                if ~(isnumeric(val) && ismatrix(val))
                    error('SSDOT:InvalidType', ...
                          'G.%s must be a numeric 2-D matrix', fname);
                end

                if isempty(val)
                    warning('SSDOT:EmptyField', 'G.%s is empty', fname);
                end
            end
        end

        function checkT(~, T)
            if ~isstruct(T)
                error('SSDOT:InvalidType', 'T must be a struct');
            end

            % Required fields
            requiredFields = ["T_HbO_brain", "T_HbR_brain", "t_HbO_brain", ...
                             "t_HbR_brain", "T_HbO_scalp", "T_HbR_scalp"];
            
            missingFields = requiredFields(~isfield(T, requiredFields));
            if ~isempty(missingFields)
                error('SSDOT:MissingFields', ...
                      'T missing required fields: %s', ...
                      strjoin(missingFields, ', '));
            end

            % Check each field
            for i = 1:numel(requiredFields)
                fieldName = requiredFields(i);
                val = T.(fieldName);

                if ~isnumeric(val)
                    error('SSDOT:InvalidType', 'T.%s must be numeric', fieldName);
                end

                if isempty(val)
                    warning('SSDOT:EmptyField', 'T.%s is empty', fieldName);
                end
            end
        end

        function checkY(~, Y)
            if ~isnumeric(Y)
                error('SSDOT:InvalidType', 'Y must be numeric');
            end

            if ~isvector(Y)
                error('SSDOT:InvalidType', 'Y must be a vector');
            end

            if isempty(Y)
                warning('SSDOT:EmptyField', 'Y is empty');
                return;
            end

            % Optional sanity checks
            if any(~isfinite(Y))
                warning('SSDOT:InvalidData', 'Y contains NaN or Inf values');
            end

            if all(Y == 0)
                warning('SSDOT:SuspiciousData', 'Y appears to be all zeros');
            end
        end

        function checkOD_SS(~, OD_SS)
            if ~(isnumeric(OD_SS) && ismatrix(OD_SS))
                error('SSDOT:InvalidType', 'OD_SS must be a numeric 2-D matrix');
            end

            if isempty(OD_SS)
                warning('SSDOT:EmptyField', 'OD_SS is empty');
            end
        end

        function checkOD_drift(~, OD_drift)
            if ~(isnumeric(OD_drift) && ismatrix(OD_drift))
                error('SSDOT:InvalidType', 'OD_drift must be a numeric 2-D matrix');
            end

            if isempty(OD_drift)
                warning('SSDOT:EmptyField', 'OD_drift is empty');
            end
        end
    end
end