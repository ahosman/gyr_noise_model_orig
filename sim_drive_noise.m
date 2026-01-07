%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   sim_drive_noise.m - Drive Channel Noise Model
%
%   Comprehensive noise model for the drive channel of an open-loop 
%   capacitive MEMS gyroscope including:
%     - Drive C/V Converter (OTA, MEMS resistor, ASIC parallel resistor)
%     - Phase-Locked Loop (PLL) components
%     - Phase Integrator
%     - Combined phase noise analysis
%
%   References: noise_model.md Section 2 (Drive Channel Noise)
%   
%   Usage: drv_noise = obj.sim_drive_noise()
%
%   Output Structure (drv_noise):
%     .cv.noise          - Drive C/V converter noise components
%     .cv.noise.vn_*     - Individual noise sources [V/√Hz]
%     .pll.noise         - PLL-related noise sources
%     .pll.num/den       - PLL transfer function coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drv_noise = sim_drive_noise(obj)
    % Initialize output structure
    drv_noise = struct();
    
    % ======================================================================
    % SECTION 1: Get Configuration and Parameters
    % ======================================================================
    
    % Get drive channel configuration
    drive = obj.imu.as.gyr.config.inf.drive;
    
    % Import physics constants
    import phys_pkg.*;
    
    % Initialize noise tracking structure
    obj.imu.as.gyr.config.inf.drive.noise = struct();
    
    % ======================================================================
    % SECTION 2: DRIVE C/V CONVERTER NOISE
    % ======================================================================
    % Section 2.2 in noise_model.md
    
    % Calculate active drive capacitance (average of DP and DN paths)
    C0_D = (obj.imu.gyro.gyro.stat.mc.Cdp + obj.imu.gyro.gyro.stat.mc.Cdn) / 2.0;
    
    % ----------------------------------------------------------------------
    % 2.1 Drive C/V OTA Noise (Section 2.2.1)
    % ----------------------------------------------------------------------
    % Formula: e_{n,out,cvd} = e_{n,opa} × (C0_D + Cp_MEMS + Cp_ASIC + Ccvd) / Ccvd
    %
    % Physical meaning: OTA input-referred noise is referred to output through
    % the impedance formed by all capacitances in the signal path.
    
    drv_noise.cv.noise.vn_cvd = ...
        obj.ctrl.noise_sel('CVD') * drive.cv.noise.Vn_opa * ...
        (C0_D + ...
         (obj.imu.gyro.gyro.par.Cpar_DP + obj.imu.gyro.gyro.par.Cpar_DN)/2 + ...
         obj.imu.as.gyr.par.Cpar_C + ...
         drive.cv.gyr_cap_cv_d) / drive.cv.gyr_cap_cv_d;
    
    % ----------------------------------------------------------------------
    % 2.2 MEMS Series Resistor Noise (Section 2.2.2)
    % ----------------------------------------------------------------------
    % Formula: e_{n,out,rs} = sqrt(2 × 4kT × Rs_D) × (C0_D + Cp_MEMS) / Ccvd
    %
    % Physical meaning: Thermal noise from electrode traces, with noise gain
    % determined only by capacitance before the resistor (feedback cap cancels).
    
    % Thermal noise density of MEMS series resistors (DP and DN combined)
    drv_noise.cv.noise.vn_rsd = ...
        obj.ctrl.noise_sel('R_MEMS_D') * ...
        sqrt(4 * phys.kT * (obj.imu.gyro.gyro.par.Rs_DP + obj.imu.gyro.gyro.par.Rs_DN));
    
    % Noise gain for MEMS resistor (only C before the Rs counts)
    drv_noise.cv.noise.rs_gain = ...
        (C0_D + (obj.imu.gyro.gyro.par.Cpar_DP + obj.imu.gyro.gyro.par.Cpar_DN)/2) / ...
        drive.cv.gyr_cap_cv_d;
    
    % Total noise density at C/V output
    drv_noise.cv.noise.vn_rs = ...
        drv_noise.cv.noise.rs_gain * drv_noise.cv.noise.vn_rsd;
    
    % ----------------------------------------------------------------------
    % 2.3 ASIC Parallel Input Resistor Noise (Section 2.2.3)
    % ----------------------------------------------------------------------
    % Formula: e_{n,out,rp} = sqrt(2 × 4kT/Rp_in) × 1/(2πf_drv × Ccvd)
    %
    % Physical meaning: Current noise from large bias resistor converted to
    % voltage through impedance of C/V integration capacitance.
    
    drv_noise.cv.noise.vn_rp = ...
        obj.ctrl.noise_sel('Rpar_D') * ...
        (sqrt(2 * 4 * phys.kT / obj.imu.as.gyr.par.Rpar_in)) / ...
        (2 * pi * obj.imu.gyro.gyro.stat.mc.fd * drive.cv.gyr_cap_cv_d);
    
    % ----------------------------------------------------------------------
    % 2.4 Total Drive C/V Converter Noise (Section 2.2.4)
    % ----------------------------------------------------------------------
    % Formula: e_{n,drive,cvd} = sqrt(vn_cvd² + vn_rs² + vn_rp²)
    %
    % RMS combination of independent noise sources
    
    drv_noise.cv.noise.vn = ...
        sqrt(drv_noise.cv.noise.vn_cvd^2 + ...
             drv_noise.cv.noise.vn_rs^2 + ...
             drv_noise.cv.noise.vn_rp^2);
    
    % Generate random noise vector (thermal white noise)
    % Noise bandwidth is fs/2 when sampled at rate fs
    drv_noise.cv.noise.vn_rand = ...
        drv_noise.cv.noise.vn * sqrt(obj.imu.as.gyr.config.inf.fs/2) * ...
        randn(1, obj.ctrl.N);
    
    % ======================================================================
    % SECTION 3: PHASE-LOCKED LOOP (PLL) NOISE
    % ======================================================================
    % Section 2.3 in noise_model.md
    
    % ----------------------------------------------------------------------
    % 3.1 Phase Integrator Signal Resistor Noise (Section 2.3.1)
    % ----------------------------------------------------------------------
    % The phase integrator uses pseudo-sinusoidal demodulation with
    % weighted input resistors for different demod phases.
    %
    % Formula: e_{n,pint,rsig} = sqrt(2×4kT × [φ_hi×R12 + φ_lo×R1] / φ_per)
    
    drv_noise.pll.noise.vn_rsig = ...
        obj.ctrl.noise_sel('Pint_Rsign') * ...
        sqrt((2 * 4 * phys.kT * drive.pll.res_sig12 * drive.pll.phi_hi / (drive.pll.phi_per/4)) + ...
             (2 * 4 * phys.kT * drive.pll.res_sig * drive.pll.phi_lo / (drive.pll.phi_per/4)));
    
    % ----------------------------------------------------------------------
    % 3.2 Phase Integrator OTA Noise (Section 2.3.2)
    % ----------------------------------------------------------------------
    % OTA noise referred to phase integrator output accounting for
    % charge integration gain.
    %
    % Formula: e_{n,pint,opa} = e_{n,opa} × K_in / K_n,in
    
    drv_noise.pll.noise.vn_int = ...
        obj.ctrl.noise_sel('Pint_opa') * drive.pll.noise.Vn_opa * ...
        drive.pll.coeff.Kin / drive.pll.coeff.Kin_noise;
    
    % ----------------------------------------------------------------------
    % 3.3 Total PLL Noise at Drive Frequency
    % ----------------------------------------------------------------------
    % Combines all contributions from C/V and phase integrator
    
    drv_noise.pll.noise.vn_pll = ...
        obj.ctrl.noise_sel('phase_noise_DRV') * ...
        sqrt(drv_noise.cv.noise.vn_cvd^2 + ...
             drv_noise.cv.noise.vn_rsd^2 + ...
             drv_noise.cv.noise.vn_rp^2 + ...
             drv_noise.pll.noise.vn_rsig^2 + ...
             drv_noise.pll.noise.vn_int^2);
    
    % ----------------------------------------------------------------------
    % 3.4 Convert to Phase Noise (dBc)
    % ----------------------------------------------------------------------
    % Conversion from voltage noise at C/V output to phase noise in dBc
    % Factor of 4/√2/√2: accounts for C/V output impedance and pseudo-sin demod
    
    drv_noise.pll.noise.phase_noise = ...
        drv_noise.pll.noise.vn_pll / ...
        (drive.pll.coeff.Kin * drive.cv.vout) / 4 / sqrt(2) / sqrt(2);
    
    drv_noise.pll.noise.phase_noise_dBc = ...
        20 * log10(drv_noise.pll.noise.phase_noise);
    
    % ======================================================================
    % SECTION 4: VCO PHASE NOISE SHAPING
    % ======================================================================
    % Section 2.4 in noise_model.md
    
    % PLL transfer function (2nd order)
    % Models loop filter response for phase noise shaping
    drv_noise.pll.num = [2 -1];      % Numerator
    drv_noise.pll.den = [1 -2 1];    % Denominator
    
    % ----------------------------------------------------------------------
    % 4.1 Generate Colored Phase Noise
    % ----------------------------------------------------------------------
    % VCO exhibits 1/f² phase noise. Generate colored noise with appropriate
    % frequency dependence.
    
    [~, phi_noise] = colouredNoise(...
        (obj.ctrl.N-1) / obj.imu.as.gyr.config.inf.fs, ...      % Duration
        obj.imu.as.gyr.config.inf.fs, ...                       % Sampling frequency
        drive.pll.noise.pll_freq, ...                           % Frequency vector
        (10.^(drive.pll.noise.pll_phnoise/20)) * ...            % Phase noise (linear)
        drive.pll.noise.scale_noise / 4 / sqrt(2) * ...         % Scaling factor
        sqrt(phys.T/300), ...                                   % Temperature compensation
        12574, 0, 0, 0, false);                                 % Random seed
    
    % Scale VCO phase noise contribution
    drv_noise.pll.noise.phase_noise_coloured = ...
        obj.ctrl.noise_sel('phase_noise_VCO') * phi_noise;
    
    % ----------------------------------------------------------------------
    % 4.2 Filter Front-End Noise Through PLL Loop Filter
    % ----------------------------------------------------------------------
    % The PLL loop filter has a low-pass characteristic. Front-end noise
    % and VCO noise combine with different transfer functions.
    
    Htest = load('filt_200hz.mat');  % Load pre-designed 200 Hz LPF
    
    % Filter C/V noise through loop filter
    drv_noise.pll.noise.phase_noise_fe_v = ...
        filter(Htest.Hlp, ...
               drv_noise.pll.noise.phase_noise_coloured * ...
               sqrt(obj.imu.as.gyr.config.inf.fs/2) .* ...
               drv_noise.cv.noise.vn_rand);
    
    % ======================================================================
    % SECTION 5: COMBINED DRIVE NOISE OUTPUT
    % ======================================================================
    
    % Optional 25 kHz residual (if present in some designs)
    drv_noise.pll.peak_at_25k = 0 * 10^(-70/20) / 64;
    
    % Total phase noise: combination of FE noise, VCO noise, and residuals
    drv_noise.pll.noise.phase_noise_drive = ...
        sqrt(2) * (...
            drv_noise.pll.peak_at_25k * ...
            sin(2*pi * obj.imu.gyro.gyro.stat.drive.tuned_frequency.mean * obj.ctrl.t) + ...
            drv_noise.pll.noise.phase_noise_fe_v + ...
            drv_noise.pll.noise.phase_noise_coloured);
    
    % ======================================================================
    % SECTION 6: SPECTRAL ANALYSIS AND VISUALIZATION
    % ======================================================================

    % Ensure backward compatibility: some control structs may not contain
    % plot_spectra_d (drive spectra flag). Default to disabled (0) if
    % the field is missing to avoid runtime errors in batch mode.
    if ~isfield(obj.ctrl, 'plot_spectra_d')
        obj.ctrl.plot_spectra_d = 0;
    end

    if obj.ctrl.plot_spectra_d == 1
        % Figure 101: Phase noise components comparison
        figure(101);
        clf;
        hold on;
        pwelch(drv_noise.pll.noise.phase_noise_drive, [], [], [], obj.imu.as.gyr.config.inf.fs);
        pwelch(drv_noise.pll.noise.phase_noise_fe_v, [], [], [], obj.imu.as.gyr.config.inf.fs);
        pwelch(drv_noise.pll.noise.phase_noise_coloured, [], [], [], obj.imu.as.gyr.config.inf.fs);
        set(gca, 'xscale', 'log');
        legend('Phase Noise Drive', 'Phase Noise FE V', 'Phase Noise Coloured');
        title('Drive Channel Phase Noise Spectrum');
        xlabel('Frequency [Hz]');
        ylabel('Power [V²/Hz]');
        grid on;
        
        % Figure 102: C/V output noise modulation
        figure(102);
        pwelch(sin(2*pi*35e3*(0:obj.ctrl.N-1)/obj.imu.as.gyr.config.inf.fs + ...
                   drv_noise.pll.noise.phase_noise_fe_v), ...
               blackman(obj.ctrl.N/8, 'periodic'), [], [], obj.imu.as.gyr.config.inf.fs);
        hold on;
        pwelch(sin(2*pi*35e3*(0:obj.ctrl.N-1)/obj.imu.as.gyr.config.inf.fs) + ...
               randn(1, obj.ctrl.N)*drv_noise.pll.noise.vn_pll/drive.cv.vout * ...
               sqrt(obj.imu.as.gyr.config.inf.fs/2), ...
               blackman(obj.ctrl.N/8, 'periodic'), [], [], obj.imu.as.gyr.config.inf.fs);
        title('Drive Signal with Noise Modulation');
        xlabel('Frequency [Hz]');
        ylabel('Power [V²/Hz]');
        grid on;
        
        % Figure 103: VCO phase noise
        figure(103);
        hold on;
        pwelch(drv_noise.pll.noise.phase_noise_coloured, [], [], [], obj.imu.as.gyr.config.inf.fs);
        title('VCO Phase Noise');
        xlabel('Frequency [Hz]');
        ylabel('Power');
        grid on;
        set(gca, 'xscale', 'log');
        
        % Figure 104: Phase noise at 400 kHz (upsampled representation)
        figure(104);
        pwelch(drv_noise.pll.noise.phase_noise_drive*16, [], [], [], ...
               obj.imu.as.gyr.config.inf.fs);
        set(gca, 'xscale', 'log');
        hold on;
        pwelch(drv_noise.pll.noise.phase_noise_drive*16, ...
               blackman(obj.ctrl.N/128, 'periodic'), [], [], ...
               obj.imu.as.gyr.config.inf.fs);
        set(gca, 'xscale', 'log');
        title('Drive Phase Noise at 400 kHz');
        xlabel('Frequency [Hz]');
        ylabel('Power');
        grid on;
    end
    
end
