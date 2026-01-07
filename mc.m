classdef mc
    properties
        imu0 (1,:) imu
        info
    end
    
    methods
        function obj = mc(nsim, model, config, mc5000File, parFile, parSheet, pmode)
            stat = gyro_mems.init_gyro_stat(nsim, mc5000File);
            stat_as = gyr.init_gyr_stat(nsim);
            fprintf("\nPreparing MC Samples ... (this can take a while)\n");
            for ii = 1:nsim
                if (ii > 1)
                    obj.imu0(ii) = imu(1, [], config, mc5000File, parFile, parSheet, pmode);
                    obj.imu0(ii).gyro = obj.imu0(ii-1).gyro;
                else
                    obj.imu0(ii) = imu(1, model, config, mc5000File, parFile, parSheet, pmode);
                end
                
                obj.imu0(ii).gyro.gyro.stat = [];
                obj.imu0(ii).gyro.gyro.stat = stat(ii);
                obj.imu0(ii).as.gyr.stat = stat_as(ii);
                [obj.imu0(ii).as.gyr.config,~] = obj.imu0(ii).as.gyr.init_gyr_trim(pmode, obj.imu0(ii).as.gyr, obj.imu0(ii).gyro);
                obj.imu0(ii).as.gyr.hist = obj.imu0(ii).as.gyr.create_stat(obj.imu0(ii).as.gyr, obj.imu0(ii).gyro);
            end
            leafNodes = mc.traverseIMUArray(obj.imu0, 'gyro.gyro.stat.mc', 'gyro');
            for idx = 1:length(leafNodes)
                fprintf('\nHandling field %s', leafNodes(idx).path);
                obj.info.gyro.param.(leafNodes(idx).path).val = leafNodes(idx).val;
            end
            x_min_mems   = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.spec_min,fieldName);
            x_max_mems   = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.spec_max,fieldName);
            x_name_mems  = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.name,fieldName);
            x_mult_mems  = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.mult,fieldName);
            x_units_mems = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.units,fieldName);
            x_color_mems = @(struct, fieldName) safe_getfield(struct.gyro.gyro.stat.color,fieldName);
            fields = fieldnames(obj.imu0(1).gyro.gyro.stat.mc);
            for i = 1:length(fields)
                try
                    % Check if field exists in the first element
                    if ~isfield(obj.imu0(1).gyro.gyro.stat.spec_min, fields{i})
                        error('Field not found');
                    end
                    obj.info.gyro.param.(['gyro__' fields{i}]).spec_min = unique(arrayfun(x_min_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyro__' fields{i}]).spec_max = unique(arrayfun(x_max_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyro__' fields{i}]).name     = unique(arrayfun(x_name_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyro__' fields{i}]).mult     = unique(arrayfun(x_mult_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyro__' fields{i}]).units    = unique(arrayfun(x_units_mems,obj.imu0, repmat(string(fields{i}), 1, nsim), 'UniformOutput', false));
                    obj.info.gyro.param.(['gyro__' fields{i}]).color    = unique(arrayfun(x_color_mems,obj.imu0, repmat(string(fields{i}), 1, nsim), 'UniformOutput', false));
                catch ME
                    fprintf(' \nWarning: Field "%s" not found in stat structures. Skipping...\n', fields{i});
                end
            end
            leafNodes = mc.traverseIMUArray(obj.imu0, 'as.gyr.stat.mc', 'gyr');
            for idx = 1:length(leafNodes)
                fprintf('\nHandling field %s', leafNodes(idx).path);
                obj.info.gyro.param.(leafNodes(idx).path).val = leafNodes(idx).val;
            end
            x_min_mems   = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.spec_min,fieldName);
            x_max_mems   = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.spec_max,fieldName);
            x_name_mems  = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.name,fieldName);
            x_mult_mems  = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.mult,fieldName);
            x_units_mems = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.units,fieldName);
            x_color_mems = @(struct, fieldName) safe_getfield(struct.as.gyr.stat.color,fieldName);
            fields = fieldnames(obj.imu0(1).as.gyr.stat.mc);
            for i = 1:length(fields)
                try
                    % Check if field exists in the first element
                    if ~isfield(obj.imu0(1).as.gyr.stat.spec_min, fields{i})
                        error('Field not found');
                    end
                    obj.info.gyro.param.(['gyr__' fields{i}]).spec_min = unique(arrayfun(x_min_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyr__' fields{i}]).spec_max = unique(arrayfun(x_max_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyr__' fields{i}]).name     = unique(arrayfun(x_name_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyr__' fields{i}]).mult     = unique(arrayfun(x_mult_mems,obj.imu0, repmat(string(fields{i}), 1, nsim)));
                    obj.info.gyro.param.(['gyr__' fields{i}]).units    = unique(arrayfun(x_units_mems,obj.imu0, repmat(string(fields{i}), 1, nsim), 'UniformOutput', false));
                    obj.info.gyro.param.(['gyr__' fields{i}]).color    = unique(arrayfun(x_color_mems,obj.imu0, repmat(string(fields{i}), 1, nsim), 'UniformOutput', false));
                catch
                    fprintf('Warning: Field "%s" not found in as.gyr.stat structures. Skipping...\n', fields{i});
                end
            end
            
            fields = fieldnames(obj.imu0(1).as.gyr.hist);
            for idx = 1:length(fields)
                fprintf('\nHandling field %s', strcat('asic__', string(fields(idx))));
                obj.info.asic.param.(strcat('asic__', string(fields(idx)))) = struct;
                obj.info.asic.param.(strcat('asic__', string(fields(idx)))).name  = obj.imu0(1).as.gyr.hist.(string(fields(idx))).name;
                obj.info.asic.param.(strcat('asic__', string(fields(idx)))).mult  = obj.imu0(1).as.gyr.hist.(string(fields(idx))).mult;
                obj.info.asic.param.(strcat('asic__', string(fields(idx)))).units = obj.imu0(1).as.gyr.hist.(string(fields(idx))).units;
                obj.info.asic.param.(strcat('asic__', string(fields(idx)))).color = obj.imu0(1).as.gyr.hist.(string(fields(idx))).color;
                for ii=1:nsim
                    if ii==1
                        obj.info.asic.param.(strcat('asic__', string(fields(idx)))).val(1)  = obj.imu0(ii).as.gyr.hist.(string(fields(idx))).val;
                    else
                        obj.info.asic.param.(strcat('asic__', string(fields(idx)))).val(end+1)  = obj.imu0(ii).as.gyr.hist.(string(fields(idx))).val;
                    end
                end
            end
            fprintf('\n');
        end
        hist_plot(obj, nestedPath);
        prob_plot(obj, nestedPath);
        scatter_plots(obj, targetPath, paramPath, thresh);
    end
    
    methods (Static, Access = public)
        leafNodes = traverseSingleStruct(currentStruct, baseName);
        leafNodes = traverseIMUArray(objArray, nestedPath, baseName);
    end
end

function value = safe_getfield(s, fieldName)
if isfield(s, fieldName)
    value = getfield(s, fieldName);
else
    value = [];
end
end
