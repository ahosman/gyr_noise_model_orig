classdef acc < handle
    properties
        % MEMS Model Parameters
        Sg = [2 4 9];           % Sensitivity values in fF/g
        beta = 0.0105;          % Non-linearity coefficient
        
        % AFE Parameters
        vcv_max = 0.85;         % Maximum CV voltage
        vsat = 100.0;           % Saturation voltage
        vref = 1.4;             % Reference voltage
        
        % Steffen AFE Model
        afe                     % AFE parameter struct from init_acc_afe
        
        % MEMS Parameters
        mems_params             % MEMS parameter struct
        
        % Simulation Parameters
        sim_params              % Simulation parameter struct
        
        % Model Configuration
        range = '32g';          % Current g-range setting
        mode = 'original';      % Operating mode
    end
    
    methods
        function obj = acc(opts)
            % Constructor: Initialize all parameter structs
            if nargin < 1
                opts = struct();
            end
            
            % Initialize MEMS parameters
            obj.mems_params = obj.get_default_mems_params();
            
            % Initialize simulation parameters  
            obj.sim_params = obj.get_default_sim_params();
            
            % Initialize AFE using Steffen model
            obj.afe = init_acc_afe(opts);
            
            % Apply trimming
            obj.afe = acfe_trim(obj.afe);
            
            % Evaluate AFE
            obj.afe = acfe_eval(obj.afe);
        end
        
        function configure_range(obj, range_str)
            % Configure accelerometer for specified g-range
            obj.range = range_str;
            
            switch range_str
                case '4g'
                    obj.afe.trm_cfb = obj.afe.trm_cfb_0;
                    obj.afe.trm_lg = obj.afe.trm_lg_0;
                case '20g'
                    obj.afe.trm_cfb = obj.afe.trm_cfb_0;
                    obj.afe.trm_lg = obj.afe.trm_lg_0;
                case '32g'
                    obj.afe.trm_cfb = obj.afe.trm_cfb_0;
                    obj.afe.trm_lg = obj.afe.trm_lg_0;
                otherwise
                    warning('Unknown range: %s', range_str);
            end
            
            % Re-evaluate AFE after range change
            obj.afe = acfe_eval(obj.afe);
        end
        
        function trim_afe(obj, opts)
            % Re-trim AFE with optional parameter overrides
            if nargin > 1
                fn = fieldnames(opts);
                for k = 1:numel(fn)
                    obj.afe.(fn{k}) = opts.(fn{k});
                end
            end
            
            obj.afe = acfe_trim(obj.afe);
            obj.afe = acfe_eval(obj.afe);
        end
        
        function noise_result = simulate_noise(obj)
            % Simulate AFE noise using Steffen model
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'analysis', 'noise'));
            noise_result = sim_acc_afe_noise(obj.afe);
        end
        
        % Main analysis methods
        acc_tf(obj)
        build_acc(obj)
        eval_nonlin_acc(obj)
        
        % Simulation methods
        out = run_simulation(obj, sim_opts)
        
        % Plotting methods
        plot_transfer_function(obj)
        plot_noise_analysis(obj)
        plot_linearity(obj)
    end
    
    methods (Access = private)
        function mems_params = get_default_mems_params(~)
            mems_params = struct();
            mems_params.m = 23.9e-9;                % proof mass
            mems_params.k = 11.98;                  % spring constant
            mems_params.b = 7.71e-5;                % damping coefficient
            mems_params.C0s = 214.6e-15;            % sensor cap
            mems_params.d = 1.30e-6;                % gap
            mems_params.aoff = 0.000;               % mechanical offset
            mems_params.C1p = 0;                    % parasitic cap parallel to C1
            mems_params.C2p = 0;                    % parasitic cap parallel to C2
            mems_params.Cpar = 1400e-15;            % parasitic cap C1/C2 to ground
            mems_params.Sn_acc = 3e-6;              % brownian noise
        end
        
        function sim_params = get_default_sim_params(~)
            fsys = 6.72e6;
            fout = 768e3;
            
            sim_params = struct();
            sim_params.Tsys = 1/fsys;
            sim_params.Ts = sim_params.Tsys / 10;
            sim_params.Tout = 1/fout;
            sim_params.pdm_clk_jitter = 0;
            sim_params.Tsys_waveform = 2;
            sim_params.Tsim = 2000e-3;
            sim_params.Talign = sim_params.Tsim/2 + 10e-6;
            sim_params.nalign = 0;
            sim_params.nplace = 4;
            sim_params.Tmodechange = sim_params.Tsim/2;
            sim_params.mode1 = 0;
            sim_params.mode2 = 0;
            
            % Timing waveform parameters
            sim_params.Tsys_bias = 1/fsys;
            sim_params.Tsys_mag = 0.01/fsys;
            sim_params.Tsys_freq = 3000.0;
            sim_params.Tsys_phase = 0.0;
            sim_params.Tsys_delay = 0.0;
            
            % Pulse parameters
            sim_params.Tsys1 = 1/fsys;
            sim_params.Tsys2 = 1.02*sim_params.Tsys1;
            sim_params.Tsys_period = sim_params.Tsim;
            sim_params.fast_length = 5000;
            
            % Acceleration waveform
            sim_params.acc_waveform = 1;
            sim_params.acc = 1.0;
            sim_params.acc_bias = 1.0;
            sim_params.acc_mag = 0.0;
            sim_params.acc_freq = 1000.0;
            sim_params.acc_phase = 0.0;
            sim_params.acc_delay = 0.0;
            sim_params.acc1 = 0.0;
            sim_params.acc2 = 1.0;
            sim_params.acc_period = sim_params.Tsim;
        end
    end
    
    methods (Static, Access = private)
        C = acc_mems(Sg, g, beta, plotfig);
        [h, Vcv_32g, Vcv_8g, acc_afe_gain_corr_32g, acc_afe_gain_corr_8g, acc_cic_gain_corr, acc_trm_dgain_s0] = c2v_acc_afe(Sg, g, C, vsat, gtrim, mode, vref, vcv_max, plotfig, color);
        h = acc_dp(range, g, vref, plotfig, color, varargin);
        [outp] = sinc3(inp, N);
        [outp] = sinc3a(inp, N);
    end
end

