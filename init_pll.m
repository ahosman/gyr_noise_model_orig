function pll = init_pll(as,cfg)
    phys=phys_pkg.phys;

    %% Phase Integrator
    pll.res_sig = sum(as.drive.fe.pll.res_sig.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_fe, 6, 1))*as.sense.gyr_runit_rppoly_scale;
    pll.res_fb  = sum(as.drive.fe.pll.res_ref.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_fe, 6, 1))*as.sense.gyr_runit_rppoly_scale;

    %% Parallel of input resistors (used for high part of pseudosin)
    pll.res_sig12 = pll.res_sig*pll.res_fb/(pll.res_sig+pll.res_fb);

    %% Ratio between the 2 levels of the pseudosin demodulation
    pll.res_ratio = 12/7;
    
    pll.phi_per = 192;    % steps into a single fdrive cycle
    phi_grd = 2;      % phi_guard duration (in 1/4 fdrive cycle)

    pll.phi_lo  = 22;     % duration of lower  step of pseudosinus demod (in 1/4 fdrive cycle)
    pll.phi_hi  = 24;     % duration of higher step of pseudosinus demod (in 1/4 fdrive cycle)

    pll.coeff.int_lo  = ((-cos(2*pi*(pll.phi_lo+phi_grd)/pll.phi_per) + cos(2*pi*phi_grd/pll.phi_per))/pll.res_ratio)/(2*pi);        % integral during pll.phi_low
    pll.coeff.int_hi  = ((-cos(2*pi*(pll.phi_hi+pll.phi_lo+phi_grd)/pll.phi_per)) + cos(2*pi*(pll.phi_lo+phi_grd)/pll.phi_per))/(2*pi);           % integral during pll.phi_high
    pll.coeff.Kin     = pll.coeff.int_lo + pll.coeff.int_hi;                                                                          % signal integration coefficient (0.1318)
    pll.coeff.Kref    = (1-(8/48))/4;                                                                                         % reference integration coefficient: guard phases every 4*fdr;

    %% Noise coefficient as calculated in bai482_Gyro_FE_concept.docx to keep in consideration demod, guard phases and folding, divided by 4 to consider 1/4th of fdrive period
    pll.coeff.Kin_noise   = 0.7613/4;
    pll.coeff.Kref_noise  = 0.9129/4;
    pll.coeff.noise_demod = 0.7613/4;

    pll.gyr_pll_cap_int = sum(as.drive.fe.pll.cap_int.*dec2binvec(cfg.otp.ana.gyr_trm_cap_prog_p, 5, 0));

    %% Noise Contribution
    switch cfg.def.corner
        case 'T'
            pll.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]            
        case 'S'
            pll.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]  
        case 'F'
            pll.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]  
    end
    
    %%VCO phase noise (fitted from measurements):
    pll.noise.scale_noise = 1;      % PLL phase noise numbers are already given for the 35kHz design for 
    pll.noise.pll_freq    = [   0  200  700  2e3  4e3 5e3 6e3  10e3 18e3 cfg.inf.fs/4 50e3  100e3 cfg.inf.fs];      % Points are taken considering a spectra with asic.fs=140kHz, the low freq part (f<12kHz doesn't change with the 25/35kHz mode, the higher freq part and the notches yes)
    pll.noise.pll_phnoise = [-117 -117 -117 -102  -93 -90 -93   -98  -94         -125 -120   -132       -136];      % Points are phase noise at 25kHz [dBc]

end

