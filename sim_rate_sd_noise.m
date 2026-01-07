%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   sim_rate_sd_noise.m - Rate Channel Sigma-Delta Modulator with LPTV Effects
%
%   Simulates the MASH 1-1 sigma-delta modulator architecture for the
%   rate output. Includes:
%     - First-stage SD (3-level) with phase noise modulation
%     - Second-stage SD (2-level) for noise shaping
%     - Digital combiner with error cancellation
%     - Quantization noise analysis
%     - Transfer function calculation
%
%   The SD modulator provides second-order noise shaping of quantization
%   noise while preserving signal integrity through digital cancellation.
%
%   LPTV Effects in SD Modulator:
%     The SD1 feedback DAC switches at drive frequency (f_d ≈ 35 kHz),
%     creating an LPTV system where drive phase noise modulates the
%     feedback timing. This causes:
%       - Multiplicative noise: FB_mod = 1 + γ*[φ_noise + I_noise/I_nom]
%       - Folding of high-frequency quantization noise into signal band
%       - Coupling of flicker noise from drive PLL into rate output
%       - Correlation between phase noise and quantization noise
%
%     The input to SD1 includes all sense channel noise (thermal, flicker,
%     Brownian) which has already been processed through the pseudo-sin
%     LPTV demodulator with coefficient K_in_noise.
%
%   References: noise_model.md Section 5 (Sigma-Delta Modulator)
%
%   Usage: [spectrum_C_rate, spectrum_S_rate, spectrum_Z_rate] = obj.sim_rate_sd_noise()
%
%   Output Structure (spectrum_*_rate):
%     .sd1.rate_out       - SD1 output spectrum [dps]
%     .sd1.qnoise         - SD1 quantization noise [dps]
%     .sd2.rate_out       - SD2 output spectrum [dps]
%     .sd.rate_out        - Final output (after digital combination) [dps]
%     [plus statistics fields from plotFunction]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectrum_C_rate, spectrum_S_rate, spectrum_Z_rate] = sim_rate_sd_noise(obj)
    
    % ======================================================================
    % SECTION 0: INITIALIZATION AND CONFIGURATION
    % ======================================================================
    
    % Control parameters
    ctrl.rate_freq = 0;  % Input rate signal frequency [Hz]
    
    % Quantizer gain factor (used in transfer function calculation)
    gainQ = 14;
    
    % Initialize output structures
    spectrum_C_rate = struct();
    spectrum_S_rate = struct();
    spectrum_Z_rate = struct();
    
    % ======================================================================
    % SECTION 1: SD INPUT SIGNAL FORMATION
    % ======================================================================
    % Section 5.2.1 in noise_model.md
    
    % SD1 input includes both signal and noise from the sense chain
    % Format: Signal = cv_gain × (rate_signal + noise)
    % The noise ratein_r_dps includes all analog noise sources (thermal,
    % flicker, mechanical Brownian, CM coupling, etc.) processed through
    % the sense channel's demodulation and integration.
    
    % C-axis: Roll rotation (rate_x)
    sd.C.u = obj.imu.as.gyr.config.inf.sense.cv.C.cv_gain * ...
        (obj.sys.rate_x.Value * cos(2*pi*ctrl.rate_freq*obj.ctrl.t) + ...
         obj.info.C.sense.noise.ratein_r_dps);
    
    % S-axis: Pitch rotation (rate_y)
    sd.S.u = obj.imu.as.gyr.config.inf.sense.cv.S.cv_gain * ...
        (obj.sys.rate_y.Value * cos(2*pi*ctrl.rate_freq*obj.ctrl.t) + ...
         obj.info.S.sense.noise.ratein_r_dps);
    
    % Z-axis: Yaw rotation (rate_z)
    sd.Z.u = obj.imu.as.gyr.config.inf.sense.cv.Z.cv_gain * ...
        (obj.sys.rate_z.Value * cos(2*pi*ctrl.rate_freq*obj.ctrl.t) + ...
         obj.info.Z.sense.noise.ratein_r_dps);
    
    % Validate SD inputs before simulation
    if any(isnan(sd.C.u(:)))
        error('sim_rate_sd_noise:InvalidInput', 'SD input for channel C contains NaN. Check ratein_r_dps or cv_gain.');
    end
    if any(isnan(sd.S.u(:)))
        error('sim_rate_sd_noise:InvalidInput', 'SD input for channel S contains NaN. Check ratein_r_dps or cv_gain.');
    end
    if any(isnan(sd.Z.u(:)))
        error('sim_rate_sd_noise:InvalidInput', 'SD input for channel Z contains NaN. Check ratein_r_dps or cv_gain.');
    end
    
    % ======================================================================
    % SECTION 2: BUILD SIGMA-DELTA SYSTEM MODELS
    % ======================================================================
    % Section 5.1 in noise_model.md
    
    % Build state-space models for the MASH 1-1 architecture
    % Returns:
    %   Hsd_model_diff: Differencing filter for error extraction
    %   Hsd_model_dig:  Digital combiner H(z) = 1 - z^(-1)
    %   ABCD_model_SD1: State-space model of first-stage SD
    %   ABCD_model_SD2: State-space model of second-stage SD
    
    [sd.Hsd_model_diff, sd.Hsd_model_dig, sd.ABCD_model_SD1, sd.ABCD_model_SD2] = ...
        obj.imu.as.gyr.rate_sd(obj);
    
    % ======================================================================
    % SECTION 3: FIRST-STAGE SIGMA-DELTA MODULATOR (SD1)
    % ======================================================================
    % Section 5.2 in noise_model.md
    
    % SD1 is a 3-level quantizer with phase noise modulation of the
    % feedback. The phase noise modulation is a critical and unique feature:
    % The SD feedback DAC switches at the drive frequency, and drive phase
    % noise modulates the switching timing, creating a multiplicative noise
    % effect that folds high-frequency quantization noise into the signal
    % band.
    
    % -----------------------------------------------------------------------
    % 3.1 FIRST-STAGE SD SIMULATION WITH PHASE NOISE MODULATION
    % -----------------------------------------------------------------------
    % Section 5.2.2-5.2.3 in noise_model.md
    %
    % LPTV Effects in Sigma-Delta Modulator:
    % The SD feedback DAC switches at the drive frequency (f_d ≈ 35 kHz),
    % making it another LPTV component in the system. Drive phase noise
    % modulates the switching timing, creating a multiplicative noise effect.
    %
    % The phase noise modulation factor is:
    % FB_mod(t) = 1 + γ × [φ_drive_noise(t) + I_DAC_noise(t)/I_DAC_nom]
    % where γ = 2/N_levels = 2/3 for 3-level SD
    %
    % This LPTV effect is critical because it:
    %   1. Modulates the quantization noise (not just signal)
    %   2. Folds high-frequency quantization noise into signal band
    %   3. Couples drive phase noise into rate output
    %   4. Creates correlation between flicker noise and quantization noise
    %
    % The modulation is passed as the last parameter to simulateDSM_mod()
    % and properly accounts for LPTV coupling of flicker noise through the
    % SD feedback path.
    
    [sd.C.v1, xn, xmax, sd.C.y1] = simulateDSM_mod(...
        sd.C.u, ...                                    % Input: rate signal + noise
        sd.ABCD_model_SD1, ...                        % State-space model
        obj.imu.as.gyr.config.inf.sense.rate.common.sd1_nlev, ... % 3-level quantizer
        [], ...                                        % No dithering
        obj.info.sense.noise.pll.sd_mod);             % Phase noise modulation
    
    % Validate SD1 output for channel C
    if any(isnan(sd.C.v1)) || any(isnan(sd.C.y1))
        error('sim_rate_sd_noise:SD1_NaNOutput', ...
            'simulateDSM_mod (SD1) returned NaN for channel C.\nInputs: u has %d NaN, pll.sd_mod has %d NaN', ...
            sum(isnan(sd.C.u(:))), sum(isnan(obj.info.sense.noise.pll.sd_mod(:))));
    end
    
    [sd.S.v1, xn, xmax, sd.S.y1] = simulateDSM_mod(...
        sd.S.u, ...
        sd.ABCD_model_SD1, ...
        obj.imu.as.gyr.config.inf.sense.rate.common.sd1_nlev, ...
        [], ...
        obj.info.sense.noise.pll.sd_mod);
    
    if any(isnan(sd.S.v1)) || any(isnan(sd.S.y1))
        error('sim_rate_sd_noise:SD1_NaNOutput', ...
            'simulateDSM_mod (SD1) returned NaN for channel S.\nInputs: u has %d NaN, pll.sd_mod has %d NaN', ...
            sum(isnan(sd.S.u(:))), sum(isnan(obj.info.sense.noise.pll.sd_mod(:))));
    end
    
    [sd.Z.v1, xn, xmax, sd.Z.y1] = simulateDSM_mod(...
        sd.Z.u, ...
        sd.ABCD_model_SD1, ...
        obj.imu.as.gyr.config.inf.sense.rate.common.sd1_nlev, ...
        [], ...
        obj.info.sense.noise.pll.sd_mod);
    
    if any(isnan(sd.Z.v1)) || any(isnan(sd.Z.y1))
        error('sim_rate_sd_noise:SD1_NaNOutput', ...
            'simulateDSM_mod (SD1) returned NaN for channel Z.\nInputs: u has %d NaN, pll.sd_mod has %d NaN', ...
            sum(isnan(sd.Z.u(:))), sum(isnan(obj.info.sense.noise.pll.sd_mod(:))));
    end
    
    % -----------------------------------------------------------------------
    % 3.2 NORMALIZE SD1 OUTPUT
    % -----------------------------------------------------------------------
    % The quantizer outputs ±2/0/-2 (3-level). Normalize to ±1/0/-1
    % for consistency with digital logic levels.
    
    sd.C.v1 = sd.C.v1 / 2;
    sd.S.v1 = sd.S.v1 / 2;
    sd.Z.v1 = sd.Z.v1 / 2;
    
    % ======================================================================
    % SECTION 4: SD1 OUTPUT ANALYSIS AND PLOTTING
    % ======================================================================
    
    % Plot SD1 rate output spectrum
    spectrum_C_rate.sd1.rate_out = plotFunction(...
        sd.C.v1 .* obj.imu.as.gyr.config.inf.sense.rate.C.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('C Total Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_S_rate.sd1.rate_out = plotFunction(...
        sd.S.v1 .* obj.imu.as.gyr.config.inf.sense.rate.S.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('S Total Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_Z_rate.sd1.rate_out = plotFunction(...
        sd.Z.v1 .* obj.imu.as.gyr.config.inf.sense.rate.Z.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('Z Total Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    % ======================================================================
    % SECTION 5: QUANTIZATION ERROR EXTRACTION AND ANALYSIS
    % ======================================================================
    % Section 5.2.4 in noise_model.md
    
    % The quantization error y1 from SD1 is the input to SD2.
    % To improve second-stage performance, a differencing filter extracts
    % a shaped version of the error that improves the MASH combination.
    
    % -----------------------------------------------------------------------
    % 5.1 EXTRACT SHAPED QUANTIZATION ERROR
    % -----------------------------------------------------------------------
    % Apply the differencing filter Hsd_model_diff to the SD1 output and
    % error to create the SD2 input u2
    
    sd.C.u2 = lsim(sd.Hsd_model_diff, [sd.C.v1 ; sd.C.y1], obj.ctrl.t)';
    sd.S.u2 = lsim(sd.Hsd_model_diff, [sd.S.v1 ; sd.S.y1], obj.ctrl.t)';
    sd.Z.u2 = lsim(sd.Hsd_model_diff, [sd.Z.v1 ; sd.Z.y1], obj.ctrl.t)';
    
    % -----------------------------------------------------------------------
    % 5.2 CALCULATE QUANTIZER GAIN FROM SD1
    % -----------------------------------------------------------------------
    % The gain_quantizer1 represents the relationship between the quantizer
    % output and the quantization error, used for transfer function analysis
    
    sd.C.gain_quantizer1 = mean(sd.C.v1 * obj.imu.as.gyr.cm.vref .* sd.C.y1 / gainQ) / ...
        mean((sd.C.y1 / gainQ).^2);
    
    sd.S.gain_quantizer1 = mean(sd.S.v1 * obj.imu.as.gyr.cm.vref .* sd.S.y1 / gainQ) / ...
        mean((sd.S.y1 / gainQ).^2);
    
    sd.Z.gain_quantizer1 = mean(sd.Z.v1 * obj.imu.as.gyr.cm.vref .* sd.Z.y1 / gainQ) / ...
        mean((sd.Z.y1 / gainQ).^2);
    
    % -----------------------------------------------------------------------
    % 5.3 PLOT SD1 QUANTIZATION NOISE
    % -----------------------------------------------------------------------
    
    spectrum_C_rate.sd1.qnoise = plotFunction(...
        sd.C.y1 .* obj.imu.as.gyr.config.inf.sense.rate.C.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('C Quantization Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_S_rate.sd1.qnoise = plotFunction(...
        sd.S.y1 .* obj.imu.as.gyr.config.inf.sense.rate.S.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('S Quantization Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_Z_rate.sd1.qnoise = plotFunction(...
        sd.Z.y1 .* obj.imu.as.gyr.config.inf.sense.rate.Z.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('Z Quantization Noise SD1 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    % ======================================================================
    % SECTION 6: SECOND-STAGE SIGMA-DELTA MODULATOR (SD2)
    % ======================================================================
    % Section 5.3 in noise_model.md
    
    % SD2 is a 2-level quantizer that shapes the quantization error from SD1.
    % Unlike SD1, SD2 uses an internal clock (not the drive frequency clock),
    % so it does NOT experience phase noise modulation from the drive.
    % This enables clean, first-order noise shaping of the SD1 error.
    
    % -----------------------------------------------------------------------
    % 6.1 SECOND-STAGE SD SIMULATION
    % -----------------------------------------------------------------------
    % Section 5.3.2 in noise_model.md
    %
    % No phase noise modulation is applied (4th parameter = [])
    
    [sd.C.v2, xn2, xmax2, y2] = simulateDSM(...
        sd.C.u2, ...                    % Input: shaped SD1 error
        sd.ABCD_model_SD2, ...          % State-space model
        2, ...                          % 2-level quantizer
        []);                            % No phase noise modulation
    
    [sd.S.v2, xn2, xmax2, y2] = simulateDSM(...
        sd.S.u2, ...
        sd.ABCD_model_SD2, ...
        2, ...
        []);
    
    [sd.Z.v2, xn2, xmax2, y2] = simulateDSM(...
        sd.Z.u2, ...
        sd.ABCD_model_SD2, ...
        2, ...
        []);
    
    % -----------------------------------------------------------------------
    % 6.2 PLOT SD2 OUTPUT SPECTRUM
    % -----------------------------------------------------------------------
    
    spectrum_C_rate.sd2.rate_out = plotFunction(...
        sd.C.v2 .* obj.imu.as.gyr.config.inf.sense.rate.C.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('C Total Noise SD2 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_S_rate.sd2.rate_out = plotFunction(...
        sd.S.v2 .* obj.imu.as.gyr.config.inf.sense.rate.S.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('S Total Noise SD2 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    spectrum_Z_rate.sd2.rate_out = plotFunction(...
        sd.Z.v2 .* obj.imu.as.gyr.config.inf.sense.rate.Z.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('Z Total Noise SD2 Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'b');
    end
    
    % ======================================================================
    % SECTION 7: DIGITAL COMBINER AND FINAL OUTPUT
    % ======================================================================
    % Section 5.4 in noise_model.md
    
    % The digital combiner applies the error cancellation filter:
    % H_dig(z) = 1 - z^(-1)
    %
    % This differentiator:
    % 1. Cancels SD1 quantization noise (via MASH architecture)
    % 2. Passes SD2 output with first-order shaping
    % 3. Results in second-order overall noise shaping: NTF = (1 - z^-1)²
    
    % -----------------------------------------------------------------------
    % 7.1 APPLY DIGITAL COMBINER FILTER
    % -----------------------------------------------------------------------
    % Section 5.4.1 in noise_model.md
    
    sd.C.v = lsim(sd.Hsd_model_dig, [sd.C.v1 ; sd.C.v2], obj.ctrl.t)';
    sd.S.v = lsim(sd.Hsd_model_dig, [sd.S.v1 ; sd.S.v2], obj.ctrl.t)';
    sd.Z.v = lsim(sd.Hsd_model_dig, [sd.Z.v1 ; sd.Z.v2], obj.ctrl.t)';
    
    % -----------------------------------------------------------------------
    % 7.2 PLOT FINAL SD OUTPUT (AFTER DIGITAL COMBINATION)
    % -----------------------------------------------------------------------
    
    spectrum_C_rate.sd.rate_out = plotFunction(...
        sd.C.v .* obj.imu.as.gyr.config.inf.sense.rate.C.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('C Total Noise SD Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'r');
    end
    
    spectrum_S_rate.sd.rate_out = plotFunction(...
        sd.S.v .* obj.imu.as.gyr.config.inf.sense.rate.S.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('S Total Noise SD Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'r');
    end
    
    spectrum_Z_rate.sd.rate_out = plotFunction(...
        sd.Z.v .* obj.imu.as.gyr.config.inf.sense.rate.Z.FS, ...
        'fs', obj.imu.as.gyr.config.inf.fs, ...
        'f0', 0, ...
        'fb', obj.ctrl.fb, ...
        'fsig', ctrl.rate_freq, ...
        'stats', true, ...
        'unit', 'dps', ...
    'plot_fft', false, ...
    'new_fig', false);
    
    if obj.ctrl.plot_fft
        set(gca, 'xscale', 'log');
        title('Z Total Noise SD Rate Out [dps]');
        ax = gca;
        h = ax.Children(4);
        set(h, 'Color', 'r');
    end
    
    % ======================================================================
    % SECTION 8: TRANSFER FUNCTION CALCULATION
    % ======================================================================
    % Section 5.4.2 in noise_model.md
    
    % Calculate the Noise Transfer Function (NTF) and Signal Transfer
    % Function (STF) for the first-stage sigma-delta modulator.
    %
    % NTF(z) = (1 - z^-1)² (second-order high-pass)
    % STF(z) = z^-1 (one sample delay)
    
    [sd.C.ntf, sd.C.stf] = calculateTF(sd.ABCD_model_SD1, sd.C.gain_quantizer1);
    [sd.S.ntf, sd.S.stf] = calculateTF(sd.ABCD_model_SD1, sd.S.gain_quantizer1);
    [sd.Z.ntf, sd.Z.stf] = calculateTF(sd.ABCD_model_SD1, sd.Z.gain_quantizer1);
    
end
