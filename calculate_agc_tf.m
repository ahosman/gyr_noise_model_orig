function [BWssensor, BW_Hz, Vacpp, xd_target, PHMARGIN_deg] = calculate_agc_tf(obj, plotfig)
    format compact
    fd     = obj.gyro.gyro.stat.mc.fd;
    Qd     = obj.gyro.gyro.stat.mc.Qd;
     
    md     = obj.gyro.gyro.stat.mc.md;
    xd     = obj.gyro.gyro.stat.mc.xd;
     
    vref   = obj.as.gyr.cm.vref;
    vcm    = obj.as.gyr.config.inf.ref.cm.vout;
    vddg   = obj.as.gyr.cm.avdd;
    vcmeff = (vcm - vddg); % Effeftive DC VCM - GVDD for the 3V drive
    
    Sda = (obj.gyro.gyro.stat.mc.Sda+obj.gyro.gyro.stat.mc.Sdi);    
    Saa = obj.gyro.gyro.stat.mc.Saa;
    
    BWssensor = fd/Qd;
    
    dxdV     = Qd*Saa*vcmeff/(pi*pi*pi*fd*fd*md);           % Drive actuation sensitivity [m/V] calculated from MEMS/ASIC parameter
    Vacpp    = xd / dxdV;                                   % Vacpp needed on AA/AI to drive MEMS at xd amplitude 
                                                            % -> MAX spec is 2.8V steady state, 3.1V startup
    gyr_cap_cv_d  = obj.as.gyr.config.inf.drive.cv.gyr_cap_cv_d; % Drive C/V capacitor [F]
    res_sig = obj.as.gyr.config.inf.drive.agc.res_sig;     % Amplitude integrator signal resisor [Ohm]
    res_ref = obj.as.gyr.config.inf.drive.agc.res_ref;     % Amplitude integrator reference resistor [Ohm]
    Cinta   = 4.807e-12;    % Amplitude integrator cap [F] 
    
    gyro_rect_integral = obj.as.gyr.config.inf.drive.agc.gyro_rect_integral;
    gyro_sine_integral = obj.as.gyr.config.inf.drive.agc.gyro_sine_integral;    

    Ksns = Sda/gyr_cap_cv_d* vcm/res_sig*gyro_sine_integral;
    Kint = 0.0;
    Kref = vref/res_ref*gyro_rect_integral;
    Kact = 2.0/(md*2*pi*fd)*Saa*vcmeff*vref/res_ref;

    xd_target = Kref*Kact/(Ksns*Kact+Kint);

    % xd_target = vref*res_sig/res_ref*gyro_rect_integral/gyro_sine_integral*gyr_cap_cv_d/(vcm*Sda);
    
    %% Transfer Functions
    
    % MEMS sensor transfer function (demodulated around fd)
    num_mems = vcm*Sda/gyr_cap_cv_d*dxdV;
    den_mems = [2*Qd/(2*pi*fd) 1];
    Hmems_tf = tf(num_mems, den_mems);
    
    Rp = 1E+11;
    num_int  = Rp/res_sig*gyro_sine_integral;
    den_int  = [Rp*Cinta 1];
    Hint_tf  = tf(num_int, den_int);
    Hint_tfz = c2d(Hint_tf, 1/fd);
    
    % Amplitude BE transfer function
    num_amp     = -2*0.6*[-10.5 10];   % typ 2x for 3v drive
    den_amp     = [1 -0.8];
    Hamp_tfz    = tf(num_amp, den_amp, 1/fd);
    Hamp_tf     = d2c(Hamp_tfz);
    
    % Drive loop total transfer function
    Htot_tf     = Hmems_tf.*Hint_tf.*Hamp_tf;
    Htot_tfz    = c2d(Htot_tf,1/fd);
    
    %% FdT Plot
    opts           = bodeoptions('cstprefs');
    opts.FreqUnits = ('Hz');
    
    if (plotfig == 1)
        figure(1)
        hold on
        bode(Htot_tf, opts)
        grid on
        
        figure(2)
        hold on
        bode(Hamp_tfz, opts)
        grid on
        
        figure(3)
        hold on
        bode(Hmems_tf, 'm', opts)
        bode(Hamp_tfz, 'r', opts)
        bode(Hint_tf, 'g', opts)
        bode(Htot_tf, 'b', opts)
        legend('MEMS + CV', 'BE Filter', 'Integrator', 'TOT_s')
        grid on
        
        figure(4)
        rlocus(Htot_tfz)    
    end

    Hloop  = feedback(Htot_tf,1);
    Hloopz = feedback(Htot_tfz,1);
    
    if (plotfig == 1)
        figure(5)
        hold on
        step(Hloopz)
    end
        
    %% Calculations
    
    W=logspace(-1, 5, 2048); % freq vector
    [MAG,PHASE] = bode(Htot_tf, W, 'b');
    imax = length(W);
    for i=1:1:imax
        magdB(i) = 20*log10(MAG(1,1,i));
        phdeg(i) = PHASE(1,1,i);
    end
    
    %%% Zero dB crossing margin calculation
    zerodB_index = min(find(magdB<0,1));
    
    % Printouts
    BW_Hz = W(zerodB_index)/(2*pi);
    PHMARGIN_deg = 180+phdeg(zerodB_index);
end
