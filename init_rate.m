function rate = init_rate(as,cfg)
    phys=phys_pkg.phys;

    rate.common.res_sig1   = sum(as.sense.fe.rate.res_sig1.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_fe, 6, 1))*as.sense.gyr_runit_rppoly_scale;    
    rate.common.res_sig2   = sum(as.sense.fe.rate.res_sig2.*dec2binvec(cfg.otp.ana.gyr_trm_res_prog_fe, 6, 1))*as.sense.gyr_runit_rppoly_scale;
    
    %% Parallel of input resistors (used for high part of pseudosin)
    rate.common.res_sig12 = rate.common.res_sig1*rate.common.res_sig2/(rate.common.res_sig1+rate.common.res_sig2);
    
    %% Ratio between the 2 levels of the pseudosin demodulation
    rate.common.res_ratio = 12/7;
    
    rate.phi_per = 192;    % steps into a single fdrive cycle
    phi_grd = 2;      % phi_guard duration (in 1/4 fdrive cycle)
    rate.phi_lo  = 22;     % duration of lower  step of pseudosinus demod (in 1/4 fdrive cycle)
    rate.phi_hi  = 24;     % duration of higher step of pseudosinus demod (in 1/4 fdrive cycle)
    
    int_lo  = ((-cos(2*pi*(rate.phi_lo+phi_grd)/rate.phi_per) + cos(2*pi*phi_grd/rate.phi_per))/rate.common.res_ratio)/(2*pi); % integral during phi_lo
    int_hi  = ((-cos(2*pi*(rate.phi_hi+rate.phi_lo+phi_grd)/rate.phi_per)) + cos(2*pi*(rate.phi_lo+phi_grd)/rate.phi_per))/(2*pi);       % integral during phi_hi
    
    rate.coeff.Kin     = int_lo + int_hi;                                                                   % signal integration coefficient (0.1318)
    rate.coeff.Kref    = (1-(8/48))/4;                                                                      % reference integration coefficient: guard phases every 4*fdr;
    
    %% Noise coefficient as calculated in bai482_Gyro_FE_concept.docx to keep in consideration demod, guard phases and folding, divided by 4 to consider 1/4th of fdrive period
    rate.coeff.Kin_noise   = 0.7613/4;
    rate.coeff.Kref_noise  = 0.9129/4;
    rate.coeff.noise_demod = 0.7613/4;

    %% Resistor used in the DAC to generate the reference current from the bandgap: 

    rate.common.gyr_cap_int = sum(as.sense.fe.rate.cap_int.*dec2binvec(cfg.otp.ana.gyr_trm_cap_prog_r,5,1));    
    rate.common.gain_70kHz  = 45;          % Integrator amplifier gain at 2*fdr [dB]
    rate.common.sd1_nlev    = 3;
    %% Rate Sigma-delta coefficients:
    
    % Coeff. 2nd order main SD
    rate.coeff.a0    = 4*rate.coeff.Kin /(cfg.inf.fs*rate.common.res_sig12*rate.common.gyr_cap_int);   
    % the 1/2 takes into account for the +2/0/-2 output of simulateDSM
    rate.coeff.f1    = 4*rate.coeff.Kref * cfg.inf.ref.idac.iout /(cfg.inf.fs*rate.common.gyr_cap_int); 
    rate.coeff.a1    = 6/6;
    % the 1/2 takes into account for the +2/0/-2 output of simulateDSM
    rate.coeff.b1    = 2/6 * cfg.inf.ref.bg.vout/2;   
    rate.coeff.gainQ = 4.5;
    
    % Coeff. 2nd order qnoise SD
    rate.coeff.c1 = 1/8 * cfg.inf.ref.bg.vout; 
    rate.coeff.c2 = 1/2;
    
    rate.coeff.f2 = 2/8 * cfg.inf.ref.bg.vout; 
    rate.coeff.a2 = 2/4;
    rate.coeff.b2 = 1/6 * cfg.inf.ref.bg.vout;
    rate.coeff.gainQ2 = 4.3348;
    
    % Coeff digital part
    rate.coeff.d2 = 2.5;
    rate.coeff.s1 = 0; 

      % Sigma delta FS
      % Rate SD Full Scale [dps]. A factor of x2 exists from original
      % equation
    rate.C.FS = (rate.coeff.Kref/rate.coeff.Kin)*2*cfg.inf.ref.idac.iout*rate.common.res_sig12/cfg.inf.sense.cv.C.cv_gain;  
    rate.S.FS = (rate.coeff.Kref/rate.coeff.Kin)*2*cfg.inf.ref.idac.iout*rate.common.res_sig12/cfg.inf.sense.cv.S.cv_gain;  
    rate.Z.FS = (rate.coeff.Kref/rate.coeff.Kin)*2*cfg.inf.ref.idac.iout*rate.common.res_sig12/cfg.inf.sense.cv.Z.cv_gain;

    %% Noise Contribution
    switch cfg.def.corner
        case 'T'
            rate.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]     -> Spot noise at 1Hz                  
        case 'S'
            rate.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]                        
        case 'F'
            rate.noise.Vn_opa = 90e-9*sqrt(phys.T/300);       % Rate int amplifier noise [V/rt(Hz)]                      
    end
end
