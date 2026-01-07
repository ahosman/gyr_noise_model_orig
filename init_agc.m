function agc = init_agc(as,cfg)    
    agc.res_sig = sum(as.drive.fe.agc.res_sig.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_fe, 6, 1));
    agc.res_ref = sum(as.drive.fe.agc.res_ref.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_a, 6, 1));

    agc.gyro_rect_integral     = 40/48;
    agc.gyro_sine_integral     = 2*(cos(4*pi/96)-cos(44*pi/96))/pi;

    agc.int_gain = (agc.gyro_rect_integral/agc.gyro_sine_integral)*...
        (agc.res_sig/agc.res_ref);

    agc.out = agc.int_gain*cfg.inf.ref.bg.vout;
end
