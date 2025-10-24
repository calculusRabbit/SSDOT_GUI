classdef ssdotClass < handle

    properties (Access = public)
        AtlasState
        Sensitivity_Matrix
        M
        A = struct( ...
            'Amatrix_brain', struct( ...
            'w1_HbO',[], ...
            'w1_HbR',[], ...
            'w2_HbO',[], ...
            'w2_HbR',[]), ...
            'Amatrix_scalp', struct( ...
            'w1_HbO',[], ...
            'w1_HbR',[], ...
            'w2_HbO',[], ...
            'w2_HbR',[]))

        T = struct( ...
            'T_HbO_brain',[],'T_HbR_brain',[], 't_HbO_brain',[],'t_HbR_brain',[], ...
            'T_HbO_scalp',[],'T_HbR_scalp',[], 't_HbO_scalp',[],'t_HbR_scalp',[])

        G = struct('G_brain',[],'G_scalp',[])

        % Parameters/users' inputs
        cfg = struct()

        % Derived vectors/matrices for modeling
        Y = [] % stacked OD vector
        OD_SS = [] % matrix
        OD_drift = [] % matrix
        H = struct('H_brain', [], 'H_scalp', [])
        b

        % container for all runs
        RunList RunClass = RunClass.empty

        % selected run
        selectedRun RunClass = RunClass.empty

        % Optional metadata
        meta struct = struct('created', datetime, 'notes', "")
    end

    methods (Access = public)
        function obj = ssdotClass()
        end

        % set value function
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
                    
                case "AtlasState"
                    obj.AtlasState = value;
                    
                case "Sensitivity_Matrix"
                    obj.Sensitivity_Matrix = value;

                case "M"
                    obj.M = value;

                case "A"
                    obj.A = value;

                case "T"           
                    obj.T = value;

                case "G"
                    obj.G = value;

                case "cfg"
                    obj.cfg = value;

                case "Y"
                    obj.Y = value;

                case "OD_SS"
                    obj.OD_SS = value;

                case "OD_drift"
                    obj.OD_drift = value;

                case "meta"
                    obj.meta = value;s
                otherwise, error("Unknown var: %s", name);
            end
        end


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
end


