%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MEMS drive TF limits for BAI430 AGC drive loop (Herschel ASIC)
%   Author: Frank Drautz
%   BAI482  First version
%   BAI430  implement as a method of class gyr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDrv_check, a_out, gain_margin] = check_drv_loop_stb(obj, plotfig)
    
    gyro = obj.gyro.gyro;
    as   = obj.as.gyr;

    m1     = gyro.mc_results.drive_mode;
    m2     = gyro.mc_results.drive.parasitic_mode.mode;
    fDrv0  = [gyro.mc_results.working_point.tuned_frequencies_with_drive(m1) gyro.mc_results.working_point.tuned_frequencies_with_drive(m2)]; % [frequncy of the drive _ frquency of the main parasitic mode]
    QDrv   = [gyro.mc_results.drive.drive_quality_factor gyro.mc_results.detection.quality_factors_all(m2)]; 
    dCadq1 = gyro.sensor_element.components.a1.capacitance_derivative(0, m1) - gyro.sensor_element.components.a2.capacitance_derivative(0, m1); % [F/(m*sqrt(kg))] (modal derivative of drive differential actuation capacitance)
    dCadq2 = gyro.sensor_element.components.a1.capacitance_derivative(0, m2) - gyro.sensor_element.components.a2.capacitance_derivative(0, m2); % [F/(m*sqrt(kg))] (modal derivative of drive differential actuation capacitance)
    dCadq  = [dCadq1 dCadq2]; % " ( " actuation " )
    
    dCddq1 = gyro.sensor_element.components.d1.capacitance_derivative(0, m1) - gyro.sensor_element.components.d2.capacitance_derivative(0, m1); % [F/(m*sqrt(kg))] (modal derivative of drive differential detection capacitance)
    dCddq2 = gyro.sensor_element.components.d1.capacitance_derivative(0, m2) - gyro.sensor_element.components.d2.capacitance_derivative(0, m2); % [F/(m*sqrt(kg))] (modal derivative of drive differential detection capacitance)
    dCddq = [ dCddq1 dCddq2 ]; % " ( " detection " )
    %% to check parameter setting
    modal_drive_displacement_for_1_m = gyro.mc_results.drive.modal_drive_displacement_for_1_m; % [m*sqrt(kg)/m] % FP.drive.modal_drive_displacement_for_1_m_with_sign
        
    %% Analysis setting
    % plotfig = 0;    % plots 0 - none (beside summary), 1 - stability related, 2 - further analysis plots
    gainVar = 1:0.05:40;
    
    paraVar = 0.88:0.02:1.2; % here multiplication of parasitic drive mode
    fDrv    = fDrv0; % reference
    % drive capacitance swing C / V actuation (resonant excitation)
    % dCddV0  = dCddq.*(0.5*dCadq*(2*(as.config.inf.ref.cm.vout-as.cm.avdd)*4/pi)).*(QDrv./(2*pi*fDrv).^2);
       
    % gain C2V -> ideal Cfb to get VBG swing (neglecting trim steps)
    % idealCfb = ( dCddq(1)*gyro.drive.amplitude.driveMean*modal_drive_displacement_for_1_m*as.config.inf.ref.cm.vout)/as.config.inf.ref.bg.vout ;
    % parasitic drive mode resonant swing V / V actuation    
    % paraDrvPeak_dVdV = (dCddV0(2)/idealCfb)*(as.config.inf.ref.cm.vout);
    
    gyro_sine_integral     = 2*(cos(4*pi/96)-cos(44*pi/96))/pi;
    
    cLegend = {};
    index_par = find(paraVar == 1);

    %find the peak closer to fd that will be useful to delete that tone and
    %calculate the gain margin of the drive
    array_parasitic = paraVar.*fDrv0(2);
    Closest_Value = interp1(array_parasitic, array_parasitic, fDrv0(1), 'nearest');
    index_closest_value = find(array_parasitic == Closest_Value);
    array_parasitic(index_closest_value) = fDrv0(1);
    for jj = 1:length(array_parasitic)
        
        % here variation of parasitic drive mode frequency
        fDrv(2) = array_parasitic(jj);
        cLegend(jj) = { [ num2str(array_parasitic(jj),'%.2f') ' Hz (freq. parasitic mode)' ] };    
        %% MEMS drive TF (incl. higher/parasitics modes)
        % MEMS drive voltage calculation (resonant excitation)
        dCddV = dCddq.*(0.5*dCadq*(2*(as.config.inf.ref.cm.vout-as.cm.avdd)*4/pi)).*(QDrv./(2*pi*fDrv).^2);
        a_out = as.config.inf.drive.agc.out*as.config.inf.drive.cv.gyr_cap_cv_d./(as.config.inf.ref.cm.vout*dCddV); % necessary drive voltage for (AGC control target assuming resonant drive mode excitation; unclear whether sensible for higher modes?)
        
        % check corresponding drive displacement [m]
        xDrv_check = (0.5*dCadq(1)*(2*(as.config.inf.ref.cm.vout-as.cm.avdd).*a_out(1)*4/pi).*(QDrv(1)/(2*pi*fDrv(1))^2))/modal_drive_displacement_for_1_m;
        
        
        %%% include trim bit here
        % fprintf('*** parameter check: \n %.2fum drive amplitude with %.2fV a_out & gyr_trm_cap_prog_cv_d = 0x%s \n', xDrv_check*1e6, a_out(1), dec2hex(as.config.otp.ana.gyr_trm_cap_prog_cv_d));
        
        %%% build MEMS TF [V/V] (i.e. CV-output Vdiff / 1V a_out drive)
        for iDrv = 1:length(fDrv)
            num_mems = (as.config.inf.ref.cm.vout/as.config.inf.drive.cv.gyr_cap_cv_d)*dCddq(iDrv).*(0.5*dCadq(iDrv)*(2*(as.config.inf.ref.cm.vout-as.cm.avdd)*4/pi));
            den_mems = [ 1 2*pi*fDrv(iDrv)/QDrv(iDrv) (2*pi*fDrv(iDrv))^2 ];
            Hmems_drv(iDrv) = tf(num_mems, den_mems);
        end
        
        Hmems_drv_tf = Hmems_drv(1);
        if jj ~= index_closest_value
            for iDrv = 2:length(fDrv)
                Hmems_drv_tf = Hmems_drv_tf + Hmems_drv(iDrv);
            end
        end
        
        %%% modulate w/ drive frequency
        s0 = 1i*2*pi*fDrv(1);
        s = tf('s');
        
        [num, den] = tfdata(Hmems_drv_tf,'v');
        
        nCoeff = length(num);
        
        % generate lower and upper side bands TF
        num_LSB = num(nCoeff);
        num_USB = num(nCoeff);
        den_LSB = den(nCoeff);
        den_USB = den(nCoeff);
        
        for order = 1:(nCoeff-1)
            num_LSB = num_LSB + num(nCoeff-order)*(s-s0)^order;
            num_USB = num_USB + num(nCoeff-order)*(s+s0)^order;
            den_LSB = den_LSB + den(nCoeff-order)*(s-s0)^order;
            den_USB = den_USB + den(nCoeff-order)*(s+s0)^order;
        end
        
        %%% sum sidebands (demodulation weighting included in int_gain)
        % corresponding to sine modulation in time-domain
        Hmems_tf = (1/(2*1i))*(...
            num_LSB/den_LSB - num_USB/den_USB);
        
        
        %%% real part -> only necessary for latter margin analysis & root-locus plots
        % seems to be numerical effect (reasoning might be that MEMS TF is real and should stay real even when modulated)
        Hmems_tf_complexe = Hmems_tf;
        [num, den] = tfdata(Hmems_tf,'v');
        Hmems_tf = tf( real(num), real(den) );
        
        %%% zpk-model necessary to avoid low frequency deviations
        % (otherwise occuring when discretizing integration of sum of LSB & USB)
        Hmems_tf = zpk(Hmems_tf);
        Hmems_tf_complexe = zpk(Hmems_tf_complexe);
        
        %% Transf functions
        fs = fDrv(1); % sampling time
        
        % saturation of gain at low f with Rp parasitic model
        Rp=1e11; % was 1e11
        
        num_int = Rp/as.config.inf.drive.agc.res_sig*gyro_sine_integral;
        den_int = [(Rp*as.drive.fe.agc.cap_int(1)) 1];
        Hint_tf = tf(num_int, den_int);
    
        Hint_tfz = c2d(Hint_tf, 1/fs);
        
        % Amplitude BE transfer function
        %num_amp = -2*0.6*[-10.5 10];       % typ 2x for 3v drive
        num_amp  = -2.0*0.6*[-10.5 10]; % typ -2x from 3v drive
        num_amp_min = -2*0.58*[-10.47 10]; % TODO: need to be fitted against MC sim - 2x for 3v drive
        num_amp_max = -2*0.6*[-10.55 10];  % TODO: need to be fitted against MC sim - 2x for 3v drive
        den_amp = [1 -0.8];
        Hamp_tfz = tf(num_amp, den_amp, 1/fs);
        Hamp_tf = d2c(Hamp_tfz);
            
        % Drive loop total transfer function
        Htot_tf = Hmems_tf.*Hint_tf.*Hamp_tf;
        Htot_tfz = c2d(Htot_tf,1/fs);
        % % alternative sample continous part & combined with BE
        Hmems_int_tfz = c2d(Hmems_tf*Hint_tf,1/fs);
        Htot_tfz  = Hmems_int_tfz*Hamp_tfz;
        
        % closed loop transfer function
        Hloop = feedback(Htot_tf,1);
        Hloopz = feedback(Htot_tfz,1);
        
        
        
        %% Analysis Plots
        if (plotfig)
            opts = bodeoptions('cstprefs');
            opts.FreqUnits = ('Hz');
            opts.Grid = 'on';
            opts.XLim = [ 1e-1, 1e6 ];
        end
        
        % stability related plots
        if plotfig == 1
            opts.XLim = [ 1e-1, 1e5 ];
            
            % OL
            fig = figure(1);
            fig.Name = 'Open loop AGC TF';
            hold on
            bode(Htot_tfz,opts);
            legend(cLegend)
            
            % CL
            fig = figure(2);
            fig.Name = 'Closed loop AGC root locus plots ';
            hold on
            rlocus(Htot_tfz,gainVar); grid on
            title([ 'AGC Root Locus (gain variation up to ' sprintf('%i',round(max(gainVar))) ')'] )
            legend(cLegend)        
        end
        
        if (plotfig > 1)
            
            %%% MEMS drive TFs
            fig = figure(1);
            fig.Name = 'MEMS drive TF';
            hold on
            bode(Hmems_drv_tf,  opts)
            % bode(Hmems_drv(1),  opts)
            legend(cLegend)
            
            fig = figure(2);
            fig.Name = 'MEMS drive modes TF';
            hold on
            for iDrv = 1:length(fDrv)
                bode(Hmems_drv(iDrv), opts)
            end
            
            fig = figure(3);
            fig.Name = 'MEMS drive demodulated TF';
            bode(Hmems_tf_complexe,opts);
            legend(cLegend)
            
            
            %%% complete TFs
            fig = figure(10);
            fig.Name = 'Overview transfer functions';
            hold on
            bode(Hmems_tf, 'm', opts)
            bode(Hamp_tfz, 'r', opts)
            bode(Hint_tf, 'g', opts)
            bode(Htot_tf, 'b', opts)
            bode(Htot_tfz, 'y', opts)
            xlim([1e-1 1e6])
            sOverview = {'MEMS + CV', 'BE Filter', 'Integrator', 'TOT_s', 'TOT_z'};
            legend( sOverview )
            
            %%% margin plots discrete / continous
            fig = figure(11);
            fig.Name = 'Total AGC OL TF (discrete)';
            hold on
            bode(Htot_tfz,opts);
            legend(cLegend)
            
            fig = figure(12);
            fig.Name = 'Total AGC OL TF (continous)';
            hold on
            bode(Htot_tf,opts);
            legend(cLegend)
            
            
            fig = figure(15);
            fig.Name = 'Total AGC OL TF (continous/discrete)';
            hold on
            bode(Htot_tf,opts);
            bode(Htot_tfz,opts);
            legend('continous', 'discrete' )
            
            
            %%% root locus -> needs real TF
            fig = figure(20);
            fig.Name = 'Root locus plots (discrete/continous)';
            subplot(2,2,1)
            rlocus(Htot_tfz)
            subplot(2,2,2)
            rlocus(Htot_tfz,1)
            subplot(2,2,3)
            rlocus(Htot_tf)
            subplot(2,2,4)
            rlocus(Htot_tf,1)
            
            %%% root locus -> gain variation
            fig = figure(21);
            fig.Name = 'Root locus plots (discrete/continous)';
            k = 1:5;
            subplot(1,2,1)
            rlocus(Htot_tfz,k); grid on
            title('Root Locus - discrete nge. fb. gain 1:5')
            subplot(1,2,2)
            k = 1:1000;
            rlocus(Htot_tf,k); grid on
            title('Root Locus - continous neg. fb. gain 1:1000')
            
            
            % pole-zero plots ( also possible for complexe/non-real TF)
            P = pzoptions;
            P.FreqUnits = 'Hz';
            P.Grid = 'on';
            
            fig = figure(22);
            fig.Name = 'Pole-zero plot (discrete)';
            hold on
            pzplot(Hloopz, P);
            xlim([-1.1 1.1]); ylim([-1.1 1.1]);
            legend(cLegend)
            %
            fig = figure(23);
            fig.Name = 'Pole-zero plot (continous)';
            hold on
            pzplot(Hloop,P);
            % xlim([-Inf 10]);
            legend(cLegend)
            
            %%% Nyquist
            fig = figure(30);
            fig.Name = 'Nyquist plot (discrete)';
            hold on
            nyquist(Htot_tfz)
            grid on
            fig = figure(31);
            fig.Name = 'Nyquist plot (continous)';
            hold on
            nyquist(Htot_tf)
            grid on
            
            
            %%% step response -> also needs real TF
            fig = figure(40);
            fig.Name = 'Step response (continous/discrete)';
            hold on
            step(Hloop)
            grid on
            step(Hloopz)
            
            fig = figure(41);
            fig.Name = 'Step response (discrete)';
            hold on
            step(Hloopz)
            grid on
            legend(cLegend)
            
            fig = figure(42);
            fig.Name = 'Step response (continous)';
            hold on
            step(Hloop)
            grid on
            legend(cLegend)
            
        end
        %% stability margin calculation (on gainVar-grid)
        polesGainVar = rlocus(Htot_tfz,gainVar);
        stabIdx(jj) = find(max(abs(polesGainVar)) < 1, 1, 'last');
        if jj == index_par
            vector = gainVar(stabIdx);
            gain_margin = vector(length(vector));
        end
    end
    
    %% plot found CL gain stability margin
    % parasitic drive mode frequency
    if (plotfig >= 1)
        fig = figure();
        fig.Name = 'CL stability gain margin (para freq)';
        plot( array_parasitic/1e3, gainVar(stabIdx), '-o', 'Linewidth', 1.5)
        hold on
        xlabel('parasitic drive mode frequency [kHz]')
        xline(fDrv0(1)/1e3, ':', 'drive', 'Linewidth', 1.5)
        xline(fDrv0(1)/1e3*1.06, ':', 'Safe Zone PLL', 'Linewidth', 1.5)
        ylabel('CL stability gain margin [-]')
        title('Found close-loop stability gain margins (from grid var.)')
        grid on
    end
    save("asic_trimming.mat")
end
