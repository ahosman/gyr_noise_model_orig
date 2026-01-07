classdef gyro_mems
    properties
        gyro
    end
    
    methods
        function obj = gyro_mems(info,model,config,mc5000File, parFile, parSheet)
            % Ensure config structure has mems field
            if ~isfield(config, 'mems')
                warning('config.mems field is missing. Creating default values.');
                config.mems = struct('quad_x', 0, 'quad_y', 0, 'quad_z', 0);
            end
            
            % Ensure all required fields in config.mems exist
            if ~isfield(config.mems, 'quad_x')
                warning('config.mems.quad_x is missing. Using default value 0.');
                config.mems.quad_x = 0;
            end
            if ~isfield(config.mems, 'quad_y')
                warning('config.mems.quad_y is missing. Using default value 0.');
                config.mems.quad_y = 0;
            end
            if ~isfield(config.mems, 'quad_z')
                warning('config.mems.quad_z is missing. Using default value 0.');
                config.mems.quad_z = 0;
            end
            
            if (~isempty(model))
                obj.gyro = gyro_mems.init_gyro(model, 'modes', 1:100, ...
                    'quadrature_dps',[config.mems.quad_x, config.mems.quad_y, config.mems.quad_z], ...
                    'drive_detection_signal', info.drive.cv.detection_target, ...
                    'cm_voltage', info.ref.cm.vout, ...
                    'calculate_mc_results', 1, ...
                    'return_sensor_element', 1);
                obj.gyro.sensor_element.config.tasks.compute_fp.optional.fast_mode = false;
                obj.gyro.fp     = gyro_mems.compute_fp(obj.gyro.sensor_element);
                [   obj.gyro.simulink_data.cp3.drive_2f, ...
                    obj.gyro.simulink_data.cn3.drive_2f, ...
                    obj.gyro.simulink_data.cp3.drive_2f_power, ...
                    obj.gyro.simulink_data.cn3.drive_2f_power, ...
                    obj.gyro.simulink_data.cp3.drive_2f_delay, ...
                    obj.gyro.simulink_data.cn3.drive_2f_delay ] = gyro_mems.init_gyro_2f(obj);
            end
            obj.gyro.stat   = gyro_mems.init_gyro_stat(1, mc5000File);
            opts = detectImportOptions(parFile);
            opts.VariableNamingRule = 'preserve';
            opts.Sheet = parSheet;
            obj.gyro.parTbl = readtable(parFile, 'Sheet', parSheet);
            obj.gyro.par    = gyro_mems.init_par(obj.gyro.parTbl);
            obj.gyro.simulink_data.parasitic_capacitances.aa_cm  = obj.gyro.par.AP_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.ai_cm  = obj.gyro.par.AN_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.da_cm  = obj.gyro.par.DP_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.di_cm  = obj.gyro.par.DN_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cp1_cm = obj.gyro.par.CGP1_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cn1_cm = obj.gyro.par.CGN1_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cp2_cm = obj.gyro.par.CGP2_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cn2_cm = obj.gyro.par.CGN2_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cp3_cm = obj.gyro.par.CGP3_CGM;
            obj.gyro.simulink_data.parasitic_capacitances.cn3_cm = obj.gyro.par.CGN3_CGM;
        end
    end
    methods (Static, Access = private)
        FP   = compute_fp(varargin);
    end
    
    methods (Static, Access = public)
        [cpz_drive_2f, cnz_drive_2f, cpz_drive_2f_power, cnz_drive_2f_power, cpz_drive_2f_delay, cnz_drive_2f_delay] = init_gyro_2f(obj);
        mems = init_gyro(varargin);
        stat = init_gyro_stat(nsim, mc5000File);
        par  = init_par_mems(parTbl);
        MEMS_tf = sense_tf(obj, plotfig)
    end
end

