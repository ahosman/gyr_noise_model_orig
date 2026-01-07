function hist = create_stat(varargin)

    as = varargin{1};
    
    if length(varargin) > 1
        gyr = varargin{2};
    end

    %% VCO and Frequencies
    hist.gyro_drv_fres.val   = as.config.inf.drive.vco.gyro_drv_fres;
    hist.gyro_drv_fres.name  = 'Drive Frequency';
    hist.gyro_drv_fres.units = 'kHz';
    hist.gyro_drv_fres.mult  = 1e-3;
    hist.gyro_drv_fres.color = 'green';

    hist.fs.val              = as.config.inf.fs;
    hist.fs.name             = 'ADC Sampling Frequency';
    hist.fs.units            = 'kHz';
    hist.fs.mult             = 1e-3;
    hist.fs.color            = 'green';

    hist.freq_init.val       = as.config.inf.drive.vco.freq_init;
    hist.freq_init.name      = 'VCO Frequency (f_{VCO})';
    hist.freq_init.units     = 'MHz';
    hist.freq_init.mult      = 1e-6;
    hist.freq_init.color     = 'green';

    hist.cap_prog_vco.val    = as.config.otp.ana.gyr_trm_cap_prog_vco;
    hist.cap_prog_vco.name   = 'gyr\_trm\_cap\_prog\_vco';
    hist.cap_prog_vco.units  = '';
    hist.cap_prog_vco.mult   = 1;
    hist.cap_prog_vco.color  = 'magenta';

    %% Drive CV
    hist.cap_prog_cv_d.val   = as.config.otp.ana.gyr_trm_cap_prog_cv_d;
    hist.cap_prog_cv_d.name  = 'gyr\_trm\_cap\_prog\_cv\_d';
    hist.cap_prog_cv_d.units = '';
    hist.cap_prog_cv_d.mult  = 1;
    hist.cap_prog_cv_d.color = 'magenta';

    hist.gyr_cap_cv_d.val   = as.config.inf.drive.cv.gyr_cap_cv_d;
    hist.gyr_cap_cv_d.name  = 'Drive Capacitance C\_CVD';
    hist.gyr_cap_cv_d.units = 'pF';
    hist.gyr_cap_cv_d.mult  = 1e12;
    hist.gyr_cap_cv_d.color = 'cyan';

    if length(varargin) > 1
        hist.vout.val   = as.config.inf.drive.cv.vout;
        hist.vout.name  = 'CVD Output';
        hist.vout.units = 'V';
        hist.vout.mult  = 1;
        hist.vout.color = 'red';
    end

    if length(varargin) > 1
        %% QC
        hist.cap_prog_qc1.val    = as.config.otp.ana.gyr_trm_cap_prog_qc_ch1;
        hist.cap_prog_qc1.name   = 'gyr\_trm\_cap\_prog\_qc\_ch1';
        hist.cap_prog_qc1.units  = '';
        hist.cap_prog_qc1.mult   = 1;
        hist.cap_prog_qc1.color  = 'magenta';
        
        hist.cap_prog_qc2.val    = as.config.otp.ana.gyr_trm_cap_prog_qc_ch2;
        hist.cap_prog_qc2.name   = 'gyr\_trm\_cap\_prog\_qc\_ch2';
        hist.cap_prog_qc2.units  = '';
        hist.cap_prog_qc2.mult   = 1;
        hist.cap_prog_qc2.color  = 'magenta';
   
        hist.cap_prog_qc3.val    = as.config.otp.ana.gyr_trm_cap_prog_qc_ch3;
        hist.cap_prog_qc3.name   = 'gyr\_trm\_cap\_prog\_qc\_ch3';
        hist.cap_prog_qc3.units  = '';
        hist.cap_prog_qc3.mult   = 1;
        hist.cap_prog_qc3.color  = 'magenta';
    
        hist.gyr_cap_ch1.val     = as.config.inf.qc.C.gyr_cap_qc;
        hist.gyr_cap_ch1.name    = 'Quadrature Compensation Capacitance C_{qc,C} - C';
        hist.gyr_cap_ch1.units   = 'pF';
        hist.gyr_cap_ch1.mult    = 1e12;
        hist.gyr_cap_ch1.color   = 'cyan';
       
        hist.gyr_cap_ch2.val     = as.config.inf.qc.S.gyr_cap_qc;
        hist.gyr_cap_ch2.name    = 'Quadrature Compensation Capacitance C_{qc,S} - S';
        hist.gyr_cap_ch2.units   = 'pF';
        hist.gyr_cap_ch2.mult    = 1e12;
        hist.gyr_cap_ch2.color   = 'cyan';
        
        hist.gyr_cap_ch3.val     = as.config.inf.qc.Z.gyr_cap_qc;
        hist.gyr_cap_ch3.name    = 'Quadrature Compensation Capacitance C_{qc,Z} - Z';
        hist.gyr_cap_ch3.units   = 'pF';
        hist.gyr_cap_ch3.mult    = 1e12;
        hist.gyr_cap_ch3.color   = 'cyan';
        
        hist.quad_comp_ch1.val   = as.config.inf.qc.C.quad_comp_dps;
        hist.quad_comp_ch1.name  = 'Quadrature Compensated (Q_{comp,C}) - C';
        hist.quad_comp_ch1.units = 'dps';
        hist.quad_comp_ch1.mult  = 1;
        hist.quad_comp_ch1.color = 'blue';
        
        hist.quad_comp_ch2.val   = as.config.inf.qc.S.quad_comp_dps;
        hist.quad_comp_ch2.name  = 'Quadrature Compensated (Q_{comp,S}) - S';
        hist.quad_comp_ch2.units = 'dps';
        hist.quad_comp_ch2.mult  = 1;
        hist.quad_comp_ch2.color = 'blue';
         
        hist.quad_comp_ch3.val   = as.config.inf.qc.Z.quad_comp_dps;
        hist.quad_comp_ch3.name  = 'Quadrature Compensated (Q_{comp,Z}) - Z';
        hist.quad_comp_ch3.units = 'dps';
        hist.quad_comp_ch3.mult  = 1;
        hist.quad_comp_ch3.color = 'blue';
        
        hist.quad_res_ch1.val    = as.config.inf.qc.C.quad_res;
        hist.quad_res_ch1.name   = 'Residual Quadrature (Q_{res,C}) - C';
        hist.quad_res_ch1.units  = 'dps';
        hist.quad_res_ch1.mult   = 1;
        hist.quad_res_ch1.color  = 'blue';
        
        hist.quad_res_ch2.val    = as.config.inf.qc.S.quad_res;
        hist.quad_res_ch2.name   = 'Residual Quadrature (Q_{res,S}) - S';
        hist.quad_res_ch2.units  = 'dps';
        hist.quad_res_ch2.mult   = 1;
        hist.quad_res_ch2.color  = 'blue';
        
        hist.quad_res_ch3.val    = as.config.inf.qc.Z.quad_res;
        hist.quad_res_ch3.name   = 'Residual Quadrature (Q_{res,Z}) - Z';
        hist.quad_res_ch3.units  = 'dps';
        hist.quad_res_ch3.mult   = 1;
        hist.quad_res_ch3.color  = 'blue';
    
        %% Sense CV
        hist.cap_prog_cv_s_ch1.val   = as.config.otp.ana.gyr_trm_cap_prog_cv_ch1;
        hist.cap_prog_cv_s_ch1.name  = 'gyr\_trm\_cap\_prog\_cv\_ch1';
        hist.cap_prog_cv_s_ch1.units = '';
        hist.cap_prog_cv_s_ch1.mult  = 1;
        hist.cap_prog_cv_s_ch1.color = 'magenta';
        
        hist.cap_prog_cv_s_ch2.val   = as.config.otp.ana.gyr_trm_cap_prog_cv_ch2;
        hist.cap_prog_cv_s_ch2.name  = 'gyr\_trm\_cap\_prog\_cv\_ch2';
        hist.cap_prog_cv_s_ch2.units = '';
        hist.cap_prog_cv_s_ch2.mult  = 1;
        hist.cap_prog_cv_s_ch2.color = 'magenta';
    
        hist.cap_prog_cv_s_ch3.val   = as.config.otp.ana.gyr_trm_cap_prog_cv_ch3;
        hist.cap_prog_cv_s_ch3.name  = 'gyr\_trm\_cap\_prog\_cv\_ch3';
        hist.cap_prog_cv_s_ch3.units = '';
        hist.cap_prog_cv_s_ch3.mult  = 1;
        hist.cap_prog_cv_s_ch3.color = 'magenta';
        
        hist.gyr_cap_cv_s_ch1.val    = as.config.inf.sense.cv.C.gyr_cap_cv_s;
        hist.gyr_cap_cv_s_ch1.name   = 'Sense Capacitance C_{cvs,C} - C';
        hist.gyr_cap_cv_s_ch1.units  = 'fF';
        hist.gyr_cap_cv_s_ch1.mult   = 1e15;
        hist.gyr_cap_cv_s_ch1.color  = 'cyan';
        
        hist.gyr_cap_cv_s_ch2.val    = as.config.inf.sense.cv.S.gyr_cap_cv_s;
        hist.gyr_cap_cv_s_ch2.name   = 'Sense Capacitance C_{cvs,S} - S';
        hist.gyr_cap_cv_s_ch2.units  = 'fF';
        hist.gyr_cap_cv_s_ch2.mult   = 1e15;
        hist.gyr_cap_cv_s_ch2.color  = 'cyan';
        
        hist.gyr_cap_cv_s_ch3.val    = as.config.inf.sense.cv.Z.gyr_cap_cv_s;
        hist.gyr_cap_cv_s_ch3.name   = 'Sense Capacitance C_{cvs,Z} - Z';
        hist.gyr_cap_cv_s_ch3.units  = 'fF';
        hist.gyr_cap_cv_s_ch3.mult   = 1e15;
        hist.gyr_cap_cv_s_ch3.color  = 'cyan';

        %% Noise Parameters (from as.stat.mc)
        % Flicker noise (K_f legacy field, now holds i_1Hz in A/sqrt(Hz))
        hist.K_f.val    = as.stat.mc.K_f;
        hist.K_f.name   = as.stat.name.K_f;
        hist.K_f.units  = as.stat.units.K_f;
        hist.K_f.mult   = as.stat.mult.K_f;
        hist.K_f.color  = as.stat.color.K_f;

        % Flicker noise (i_1Hz new field name)
        hist.i_1Hz.val   = as.stat.mc.i_1Hz;
        hist.i_1Hz.name  = as.stat.name.i_1Hz;
        hist.i_1Hz.units = as.stat.units.i_1Hz;
        hist.i_1Hz.mult  = as.stat.mult.i_1Hz;
        hist.i_1Hz.color = as.stat.color.i_1Hz;

        % Thermal noise floor
        hist.v_th.val    = as.stat.mc.v_th;
        hist.v_th.name   = as.stat.name.v_th;
        hist.v_th.units  = as.stat.units.v_th;
        hist.v_th.mult   = as.stat.mult.v_th;
        hist.v_th.color  = as.stat.color.v_th;

    end

end

