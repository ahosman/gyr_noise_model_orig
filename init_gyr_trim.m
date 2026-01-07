function [as_config, gyro_config] = init_gyr_trim(pmode, varargin)

    as_config   = struct;
    gyro_config = struct;
    
    asic_struct = varargin{1};
    if length(varargin) > 1
        gyro_struct = varargin{2}.gyro;
    end

    as_config.otp = gyr.init_gyr_otp();
    
    as_config.inf.FSdps        = 6000.0;

    as_config.inf.afe_clkscheme_common_offset = 6;

    %% Definitions   
    as_config.def.corner     = 'T';
    as_config.def.mode       = pmode;
    as_config.def.s2fc       = 0;           % 2f-signal compensation in C/V input IS NOT implemented
    
    %% Reference
    [as_config.inf.ref, as_config.otp.ana.gyr_trm_res_prog_i, as_config.otp.ana.gyr_trm_res_prog_bg]  = gyr.init_ref(asic_struct,as_config);

    %% VCO
    if length(varargin) == 1
        [as_config.inf.drive.vco,as_config.otp.ana.gyr_trm_cap_prog_vco ] = gyr.init_vco(asic_struct, as_config);
    else
        [as_config.inf.drive.vco, as_config.otp.ana.gyr_trm_cap_prog_vco ] = gyr.init_vco(asic_struct, as_config, gyro_struct);
    end
    as_config.inf.fs        = 4*as_config.inf.drive.vco.gyro_drv_fres; % sample frequency for ADC

    %% AGC
    as_config.otp.ana.gyr_trm_res_prog_a  = as_config.otp.ana.gyr_trm_res_prog_i;
    as_config.otp.ana.gyr_trm_res_prog_fe = as_config.otp.ana.gyr_trm_res_prog_i;
    as_config.inf.drive.agc = gyr.init_agc(asic_struct,as_config);
    
    %% PLL
    as_config.inf.drive.pll = gyr.init_pll(asic_struct,as_config);
    
    %% Drive C2V
    if length(varargin) > 1
        [gyro_config.inf.drive.cv, ...
         as_config.inf.drive.cv, ...
         as_config.otp.ana.gyr_trm_cap_prog_cv_d ] = ...
         gyr.init_drive_cv(asic_struct, as_config, gyro_struct);
    else
        [gyro_config, as_config.inf.drive.cv, ...
         as_config.otp.ana.gyr_trm_cap_prog_cv_d ] = ...
         gyr.init_drive_cv(asic_struct, as_config);
    end

    %% Sense C2V
    if length(varargin) > 1
        [as_config.inf.sense.cv, ...
         as_config.otp.ana.gyr_trm_cap_prog_cv_ch1, ...
         as_config.otp.ana.gyr_trm_cap_prog_cv_ch2, ...
         as_config.otp.ana.gyr_trm_cap_prog_cv_ch3 ] = ...
         gyr.init_sense_cv(asic_struct, as_config, gyro_struct);
    end
    
    %% Quadrature compensation
    if length(varargin) > 1
        [as_config.inf.qc, ...
            as_config.otp.ana.gyr_trm_cap_prog_qc_ch1, ...
            as_config.otp.ana.gyr_trm_cap_prog_qc_ch2, ...
            as_config.otp.ana.gyr_trm_cap_prog_qc_ch3 ] = ...
            gyr.init_qc(asic_struct, as_config, gyro_struct);
    end
    
    %% Rate ADC Integrator
    if length(varargin) > 1
        as_config.inf.sense.rate = gyr.init_rate(asic_struct, as_config);
    %else
    end
            
    if length(varargin) > 1
        as_config.otp.gyr_afe.gyr_notch_delay_ch1 = abs(round(2*gyro_struct.stat.mc.fd./gyro_struct.stat.mc.fsplit_C));
        as_config.otp.gyr_afe.gyr_notch_delay_ch2 = abs(round(2*gyro_struct.stat.mc.fd./gyro_struct.stat.mc.fsplit_S));
        as_config.otp.gyr_afe.gyr_notch_delay_ch3 = abs(round(2*gyro_struct.stat.mc.fd./gyro_struct.stat.mc.fsplit_Z));
    end
end
