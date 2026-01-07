classdef noise_model < handle
    properties
        ctrl        % Simulation control parameters
        imu         % IMU/gyroscope configuration object
        sys         % System state (rates, temperature, etc.)
        info        % Noise simulation results
        spectrum    % Spectral analysis results
        adev        % Allan deviation results
        config      % Configuration metadata
    end
    
    methods
        % ====================================================================
        % CONSTRUCTOR
        % ====================================================================
        function obj = noise_model(imu, sys, noise_sel)
            % NOISE_MODEL Constructor - Initialize and run noise simulation
            %
            % Inputs:
            %   imu       - IMU configuration object
            %   sys       - System state structure
            %   noise_sel - Noise configuration (string or function handle)
            
            % Store inputs
            obj.imu = imu;
            obj.sys = sys;
            
            % Initialize control parameters (static method call)
            obj.ctrl = noise_model.init_ctrl();
            
            % Initialize noise configuration
            if ischar(noise_sel) || isstring(noise_sel)
                % String preset (e.g., 'All', 'Thermal_Only')
                obj.ctrl.noise_sel = noise_model.init_noise_cfg(char(noise_sel));
                obj.config.mode = char(noise_sel);
            elseif isa(noise_sel, 'function_handle')
                % Pre-configured function handle
                obj.ctrl.noise_sel = noise_sel;
                obj.config.mode = 'Custom';
            else
                error('noise_sel must be string or function handle');
            end
            
            % Store configuration metadata
            obj.config.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            obj.config.fs = obj.imu.as.gyr.config.inf.fs;
            obj.config.N = obj.ctrl.N;
            obj.config.duration = obj.ctrl.N / obj.imu.as.gyr.config.inf.fs;
            
            % Generate time vector
            obj.ctrl.t = 0:1/obj.imu.as.gyr.config.inf.fs:(obj.ctrl.N-1)/obj.imu.as.gyr.config.inf.fs;
            
            % Run noise simulations
            fprintf('\n=== Starting Noise Model Simulation ===\n');
            fprintf('Configuration: %s\n', obj.config.mode);
            fprintf('Duration: %.2f s (%d samples @ %.0f kHz)\n', ...
                    obj.config.duration, obj.config.N, obj.config.fs/1000);
            
            % Drive channel
            fprintf('Simulating drive channel noise...\n');
            obj.info.drive.noise = obj.sim_drive_noise();
            
            % Sense channel
            fprintf('Simulating sense channel noise...\n');
            [obj.info.sense.noise, ...
             obj.info.C.sense.noise, ...
             obj.info.S.sense.noise, ...
             obj.info.Z.sense.noise, ...
             obj.spectrum.C.sense, ...
             obj.spectrum.S.sense, ...
             obj.spectrum.Z.sense] = obj.sim_sense_noise();
            
            % Sigma-delta modulator
            fprintf('Simulating sigma-delta modulator...\n');
            [obj.spectrum.C.rate, ...
             obj.spectrum.S.rate, ...
             obj.spectrum.Z.rate] = obj.sim_rate_sd_noise();
            
            % Allan deviation (optional - can be time consuming)
            if obj.ctrl.calc_adev
                fprintf('Calculating Allan deviation...\n');
                [obj.adev.C, obj.adev.S, obj.adev.Z] = obj.allandev();
            end
            
            fprintf('=== Simulation Complete ===\n\n');
            
            % Print quick summary
            obj.print_summary();
        end
        
        % ====================================================================
        % PUBLIC METHODS
        % ====================================================================
        
        % Method declarations - implementations are in separate files
        noise_rms = get_total_noise(obj, axis)
        report = generate_report(obj, options)
        validate(obj)
        print_summary(obj)
        plot_spectra(obj, axis)
        export_data(obj, filename)
        
        % ====================================================================
        % SIMULATION METHODS
        % ====================================================================
        
        % These call external sim_*.m files
        drv_noise = sim_drive_noise(obj)
        [sns_noise, C_sense, S_sense, Z_sense, ...
         spectrum_C_cv, spectrum_S_cv, spectrum_Z_cv] = sim_sense_noise(obj)
        [spectrum_C_rate, spectrum_S_rate, spectrum_Z_rate] = sim_rate_sd_noise(obj)
        
        % ====================================================================
        % ANALYSIS METHODS
        % ====================================================================
        
        [adev_C, adev_S, adev_Z] = allandev(obj)
    end
    
    methods (Static)
        % ====================================================================
        % STATIC METHODS
        % ====================================================================
        
        % Initialization methods - implementations in separate files
        ctrl = init_ctrl()
        noise_sel = init_noise_cfg(noise_cfg_sel)
        
        % Utility method
        w = pseudo_sin(t, N)
    end
    
    methods (Access = private)
        % ====================================================================
        % PRIVATE HELPER METHODS
        % ====================================================================
        
        % Implementations in separate files
        adev_result = calc_single_axis_adev(obj, signal)
        print_noise_report(obj, report)
        save_report_file(obj, report, filename)
    end
end
