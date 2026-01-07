function [q_d, t_d] = calculate_residual_ampl(obj)
    gyro = obj.gyro.gyro;
    as   = obj.as.gyr;
    FP   = gyro.fp;

    % default values
    use_freq_nl = 1; % should be always used since tauVGM<<tauDrv

    % drive mode id
    idDrv = FP.drive_mode;

    % raw circ freq [rad/s]
    if isfield(FP.detection, 'raw_frequency')
        omDrv_r = (2*pi)*FP.detection.raw_frequency.drive;
    else
        omDrv_r = (2*pi)*FP.working_point.raw_frequency.drive;
    end
    if use_freq_nl
        if isfield(FP.detection,'raw_frequency_with_drive')
            omDrv_r = (2*pi) * FP.detection.raw_frequency_with_drive.drive;
        else
            omDrv_r = (2*pi) * FP.working_point.raw_frequency_with_drive.drive;
        end
    end
    % tun circ freq [rad/s]
    if use_freq_nl
        if isfield(FP.detection,'tuned_frequency_with_drive')
            omDrv_t = (2*pi) * FP.detection.tuned_frequency_with_drive.drive;
        else
            omDrv_t = (2*pi) * FP.working_point.tuned_frequency_with_drive.drive;
        end
        if isfield(FP.detection,'pll_drive_frequency')
            omDrv_t = (2*pi) * FP.drive.pll_drive_frequency;
        end
    end
    
    % quality factor
    Qfac_d = FP.detection.quality_factor.drive;
    
    % steady modal ampl and phase before shut down
    % qDrv = FP.drive.modal_displacement_with_sign;
    % pDrv = -pi/2*0;
    qDrv = abs(FP.drive.modal_drive_vector(idDrv));
    pDrv = phase(FP.drive.modal_drive_vector(idDrv));
    
    % damped circular frequency
    % Omega_d = omDrv_r * sqrt(1-1/(2*Qfac_d)^2);
    Omega_d = omDrv_r * sqrt(1-(omDrv_r/omDrv_t)^2/(2*Qfac_d)^2);
    
    % decay time constant [s]
    tau_d = 2*Qfac_d / omDrv_r;
    
    % modal ampl and vel at shut down
    tShut = rand(1)*(2*pi)/omDrv_t; % pick random instant in drive period
    qDrv_0 = qDrv * cos(omDrv_t * tShut + pDrv);
    qpDrv_0 = -omDrv_t * qDrv * sin(omDrv_t * tShut + pDrv);
    
    % modal equivalent amplitude and phase at shut down
    cosC = qDrv_0;
    sinC = -(1/Omega_d) * (qDrv_0/tau_d + qpDrv_0);
    ampR_d = sqrt(cosC^2 + sinC^2);
    phiR_d = atan2(sinC, cosC);

    % stiffness time decay constant [s] (notice that stiff goes with VGM^2  where VGM=VGM0*exp(-t/tauVGM)then 1/2 coef is needed in tauK)
    if isfield(FP,'mcConfig')
        if ~isfield(FP.mcConfig.asic,'tau_cm')
            FP.mcConfig.asic.tau_cm = 2e-3;
            warning('Using default value: FP.mcConfig.asic.tau_cm = 2e-3')
        end
        tau_k = .5 * FP.mcConfig.asic.tau_cm;
    else
        if ~isfield(FP.misc,'tau_cm')
            FP.misc.tau_cm = 2e-3;
            warning('Using default value: FP.misc.tau_cm = 2e-3')
        end
        tau_k = .5 * FP.misc.tau_cm;
    end
    nom = gyroObj.settings.simulation.number_of_modes;
    qd  = zeros(length(tDecay));
    qd  = ampR_d*exp(-tDecay/tau_d);

end

