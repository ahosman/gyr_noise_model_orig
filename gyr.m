classdef gyr < handle
    properties
        cm
        drive
        sense
        noise
        config
        par
        stat
        hist
    end
    
    methods
        function obj = gyr(pmode)
            gyr.init_gyr(obj);
            obj.stat = gyr.init_gyr_stat(1);
            [obj.config,~] = gyr.init_gyr_trim(pmode, obj);
            obj.hist = gyr.create_stat(obj);
            obj.par = gyr.init_par();
            gyr.init_gyr_noise(obj);
        end
    end
    methods (Static, Access = private)
        init_gyr(obj);
        init_gyr_noise(obj);
        par = init_par();
        otp = init_gyr_otp();
        [ref, res_prog_i, res_prog_bg] = init_ref(as,cfg);
        [vco, cap_prog_vco] = init_vco(varargin);
        drv_agc = init_agc(as,cfg);
        [drv_pll, pll_hist] = init_pll(as,cfg);
        [gyr_drv, as_drv_cv, cap_prog_cv_d] = init_drive_cv(varargin);
        [qc, cap_prog_qc1, cap_prog_qc2, cap_prog_qc3] = init_qc(varargin);
        [as_sns_cv, cap_prog_cv_s_ch1, cap_prog_cv_s_ch2, cap_prog_cv_s_ch3]  = init_sense_cv(varargin);
        rate = init_rate(as,cfg);
    end
    methods (Access = public)
    end
    methods (Static, Access = public)
        stat = init_gyr_stat(nsim);
        hist = create_stat(varargin);
        [as_config, gyr_config] = init_gyr_trim(pmode, varargin);
        H_cic = cic_tf(f,fs,R,ord,plotfig);
        fun_folding = calc_folding (f, fs, R);
        [ Hsd_model_SD1, Hsd_model_SD2, Hsd_model_diff, Hsd_model_dig,  ABCD_model_SD1, ABCD_model_SD2 ] ...
            = rate_sd(obj);
        cic_gain_corr = calc_cic_gain_corr(N, M, cic_osr_lut);
    end
end

