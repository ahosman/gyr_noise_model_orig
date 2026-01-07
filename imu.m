classdef imu
    properties
        as asic
        gyro gyro_mems
    end

    methods
        function obj = imu(nstr, model, config, mc5000File, parFile, parSheet, pmode)
            obj.as   = asic(nstr, pmode);
            obj.gyro = gyro_mems(obj.as.gyr.config.inf, model, config, ...
                mc5000File, parFile, parSheet);
            [obj.as.gyr.config,~] = obj.as.gyr.init_gyr_trim(pmode, obj.as.gyr, obj.gyro);
            obj.gyro = gyro_mems(obj.as.gyr.config.inf, model, config, ...
                mc5000File, parFile, parSheet);
            obj.as.gyr.hist = obj.as.gyr.create_stat(obj.as.gyr);
            [obj.as.gyr.config,~] = obj.as.gyr.init_gyr_trim(pmode, obj.as.gyr, obj.gyro);
        end
        [BWssensor, BW_Hz, Vacpp, xd_target, PHMARGIN_deg] = calculate_agc_tf(obj, plotfig);
        [Sdd_C, Sdd_S, Sdd_Z, stot_C, stot_S, stot_Z, gyr_trm_dgain_s0_ch1, gyr_trm_dgain_s0_ch2, gyr_trm_dgain_s0_ch3] = calculate_sens(obj)
        [xDrv_check, a_out, gain_margin] = check_drv_loop_stb(obj, plotfig);
        gyr_tf(obj, gyr_mode, filt_odr, filt_type, filt_bw, filt_rolloff_filt_notch_bypass);
    end
end
