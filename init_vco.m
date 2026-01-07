function [vco, cap_prog_vco] = init_vco(varargin)
    as = varargin{1};
    cfg = varargin{2};
    if length(varargin) > 2
        gyro = varargin{3};
    end
    vco.n_cycle_afe  = 192;
    %% Frequencies
    if length(varargin) > 2
        vco.gyro_drv_fres = gyro.stat.mc.fd;
    else
        vco.gyro_drv_fres = 35.0E+03;
    end
    vco.fosc_nom    = 6.54E+06; % nominal FOSC frequency at trim = 0
    vco.fosc_target = 6.72E+06; % target FOSC frequeny (for FOSC_freq_offset = 0.0 and gyr_drv_fres = 35kHz)

    vco.fsys_target    = vco.n_cycle_afe*vco.gyro_drv_fres;

    vco.fosc_delta     = vco.fsys_target-vco.fosc_nom;
    vco.ib_fosc_target = vco.fosc_delta*1.24E-09/21.1E+03 + cfg.inf.ref.ib.fosc_nom;

    vco.res_sig   = sum(as.drive.vco.res_sig.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_i, 6, 0));
    vco.res_ref   = sum(as.drive.vco.res_ref.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_vco, 6, 0));
    vco.cap_sig   = -0.5*(as.drive.vco.i_curr)/(-1.0*as.drive.vco.i_curr*vco.res_ref*vco.fsys_target);
    cap_prog_vco = ...
        round((vco.cap_sig-as.drive.vco.cap_vco(1)-as.drive.vco.cap_par)/as.drive.vco.cap_vco(2));
    vco.freq_init = ...
        -0.5*(as.drive.vco.i_curr)/(sum(as.drive.vco.cap_vco.*dec2binvec(cap_prog_vco, 5, 0))+as.drive.vco.cap_par)...
        *-1/(as.drive.vco.i_curr*sum(as.drive.vco.res_ref.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_vco, 6, 0)));
end
