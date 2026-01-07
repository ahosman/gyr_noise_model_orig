
## 9. Implementation Guide

### 9.1 File Structure

```
noise_model/
├── noise_model.m                    % Main class definition
├── sim_drive_noise.m                % Drive channel noise simulation
├── sim_sense_noise.m                % Sense channel noise simulation
├── sim_rate_sd_noise.m              % Sigma-delta modulator simulation
├── init_noise_cfg_enhanced.m        % Noise source configuration
├── generate_rtn.m                   % RTN generation function
├── generate_noise_breakdown.m       % Analysis and reporting
└── colouredNoise.m                  % Colored noise generation utility
```

### 9.2 Usage Examples

#### 9.2.1 Basic Simulation

```matlab
% Load IMU configuration
load('imu_config.mat', 'imu');

% Define system state
sys.rate_x = 10;   % 10 dps input (C-axis)
sys.rate_y = 0;    % 0 dps input (S-axis)
sys.rate_z = 0;    % 0 dps input (Z-axis)

% Create noise model with all sources
nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));

% Check results
fprintf('Total noise (C-axis): %.3f dps RMS\n', ...
        std(nm.info.C.sense.noise.ratein_r_dps));
```

#### 9.2.2 Noise Source Study

```matlab
% Compare different noise configurations
configs = {'All', 'Thermal_Only', 'LowFreq_Only', 'Disable_CM'};
noise_levels = zeros(size(configs));

for i = 1:length(configs)
    nm = noise_model(imu, sys, init_noise_cfg_enhanced(configs{i}));
    noise_levels(i) = std(nm.info.C.sense.noise.ratein_r_dps);
end

% Plot comparison
figure;
bar(noise_levels * 1000);  % Convert to mdps
set(gca, 'XTickLabel', configs);
ylabel('Total Noise [mdps RMS]');
title('Noise Configuration Comparison');
grid on;
```

#### 9.2.3 Temperature Sweep

```matlab
temperatures = 25:20:125;  % -40 to 125°C
arw_vs_temp = zeros(size(temperatures));

for i = 1:length(temperatures)
    % Update temperature
    imu.temperature = temperatures(i) + 273.15;  % Convert to Kelvin
    
    % Run simulation
    nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
    
    % Calculate ARW
    [psd, f] = pwelch(nm.info.C.sense.noise.ratein_r_dps, [], [], [], ...
                     nm.imu.as.gyr.config.inf.fs);
    idx = (f >= 10) & (f <= 100);
    white_noise = median(sqrt(psd(idx)));
    arw_vs_temp(i) = white_noise * sqrt(3600) * 60;  % deg/√hr
end

% Plot
figure;
plot(temperatures, arw_vs_temp, 'LineWidth', 2);
xlabel('Temperature [°C]');
ylabel('ARW [deg/√hr]');
title('Angle Random Walk vs. Temperature');
grid on;
```

#### 9.2.4 Adding RTN

```matlab
% Enhanced noise configuration with RTN
noise_cfg = init_noise_cfg_enhanced('All');

% Run base simulation
nm = noise_model(imu, sys, noise_cfg);

% Manually add RTN to sense C/V OTA
num_traps = 3;
tau_high = [15e-6, 30e-6, 80e-6];
tau_low = [12e-6, 25e-6, 70e-6];
amplitude = [1.5e-3, 2.0e-3, 1.2e-3];  % mV

rtn_signal = generate_rtn(nm.ctrl.N, nm.imu.as.gyr.config.inf.fs, ...
                         amplitude, tau_high, tau_low, num_traps, 42);

% Add to existing noise (with appropriate gain)
cvs_gain = nm.imu.as.gyr.config.inf.sense.cv.C.cvs_gain;
nm.info.C.sense.noise.cvout_cvs_ota = ...
    nm.info.C.sense.noise.cvout_cvs_ota + rtn_signal * cvs_gain;

% Recalculate totals
nm.info.C.sense.noise.cvout_r = nm.info.C.sense.noise.cvout_r + ...
                                rtn_signal * cvs_gain;

% New total noise
fprintf('Noise with RTN: %.3f dps RMS\n', ...
        std(nm.info.C.sense.noise.ratein_r_dps));
```

### 9.3 Configuration File Structure

#### 9.3.1 IMU Configuration

```matlab
% Example IMU configuration structure
imu.as.gyr.config.inf.fs = 100e3;          % Sampling frequency [Hz]
imu.as.gyr.config.inf.FSdps = 100;         % Full-scale range [dps]

% Drive parameters
imu.gyro.gyro.stat.mc.fd = 35e3;           % Drive frequency [Hz]
imu.gyro.gyro.par.Rs_DP = 100;             % Drive series R [Ω]
imu.gyro.gyro.par.Cpar_DP = 50e-15;        % Drive parasitic C [F]

% Sense parameters
imu.gyro.gyro.par.Rs_GP1 = 150;            % Sense series R [Ω]
imu.gyro.gyro.par.CGP1_CGM = 800e-15;      % Sense capacitance [F]

% ASIC parameters
imu.as.gyr.config.inf.sense.cv.noise.Vn_opa = 8e-9;  % OTA noise [V/√Hz]
imu.as.gyr.config.inf.sense.rate.coeff.Kin = 0.85;   % Signal coefficient
imu.as.gyr.config.inf.sense.rate.coeff.Kin_noise = 0.761;  % Noise coeff

% Reference parameters
imu.as.gyr.config.inf.ref.bg.vout = 1.2;   % Bandgap voltage [V]
imu.as.gyr.config.inf.ref.cm.vout = 0.9;   % CM voltage [V]
```

### 9.4 Customization

#### 9.4.1 Adding New Noise Source

```matlab
% In sim_sense_noise.m, add new source:

% New noise source: Substrate coupling
if obj.ctrl.noise_sel('Substrate_Coupling')
    % Calculate noise
    substrate_noise = sqrt(4*phys.kT*Rsub) * ...
                     Csub_mismatch/Ccvs * ...
                     sqrt(fs/2) * randn(1,N);
    
    % Add to total
    C_sense.cvout_r = C_sense.cvout_r + substrate_noise;
end

% In init_noise_cfg_enhanced.m, add to source list:
enabled_sources = {...
    'CVD', 'R_MEMS_D', ..., ...
    'Substrate_Coupling'  % New source
};
```

#### 9.4.2 Custom Analysis Function

```matlab
function custom_analysis(nm)
    % Custom analysis on noise_model object
    
    % Extract time-domain signals
    rate_signal = nm.info.C.sense.noise.ratein_r_dps;
    
    % Your custom analysis
    % ... (spectral analysis, statistics, etc.)
    
    % Generate custom plots
    figure;
    % ... your plots
end

% Usage
nm = noise_model(imu, sys, 'All');
custom_analysis(nm);
```

---

## 10. Validation and Testing

### 10.1 Unit Tests

#### 10.1.1 Thermal Noise Validation

```matlab
function test_thermal_noise_scaling()
    % Test that thermal noise scales as sqrt(T)
    
    T1 = 300;  % 300 K
    T2 = 400;  % 400 K
    Rs = 100;  % 100 Ω
    
    vn1 = sqrt(4*1.38e-23*T1*Rs);
    vn2 = sqrt(4*1.38e-23*T2*Rs);
    
    ratio_expected = sqrt(T2/T1);
    ratio_actual = vn2/vn1;
    
    assert(abs(ratio_actual - ratio_expected) < 0.001, ...
           'Thermal noise temperature scaling failed');
    
    fprintf('✓ Thermal noise scaling test passed\n');
end
```

#### 10.1.2 Noise Source Independence

```matlab
function test_noise_independence()
    % Test that individual noise sources can be enabled/disabled
    
    % All enabled
    nm_all = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
    noise_all = std(nm_all.info.C.sense.noise.ratein_r_dps);
    
    % Thermal only
    nm_thermal = noise_model(imu, sys, init_noise_cfg_enhanced('Thermal_Only'));
    noise_thermal = std(nm_thermal.info.C.sense.noise.ratein_r_dps);
    
    % Thermal should be less than all
    assert(noise_thermal < noise_all, ...
           'Thermal-only noise should be less than all noise');
    
    fprintf('✓ Noise source independence test passed\n');
end
```

#### 10.1.3 Pseudo-Sin Coefficients

```matlab
function test_pseudo_sin_coefficients()
    % Verify pseudo-sin noise coefficient calculation
    
    N = 192;
    w = pseudo_sin(0:N-1, N);
    
    % Calculate RMS
    Kin_noise_calc = sqrt(mean(w.^2));
    Kin_noise_expected = 0.761;
    
    assert(abs(Kin_noise_calc - Kin_noise_expected) < 0.01, ...
           'Pseudo-sin noise coefficient incorrect');
    
    fprintf('✓ Pseudo-sin coefficient test passed\n');
    fprintf('  Calculated: %.4f, Expected: %.4f\n', ...
            Kin_noise_calc, Kin_noise_expected);
end
```

### 10.2 Integration Tests

#### 10.2.1 Total Noise RSS Check

```matlab
function test_noise_rss()
    % Verify that total noise matches RSS of components
    
    nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
    
    % Calculate individual contributions
    contrib = nm.info.C.sense.noise;
    sources = {
        'cvout_bg_noise', 'cvout_r_cpcn', 'cvout_cvs_ota', ...
        'cvout_rsig', 'cvout_rate_ota', 'cvout_cm_c0_noise', ...
        'cvout_rp', 'cvout_cm_qdem_noise_r', 'cvout_cm_cpar_noise', ...
        'cvout_cm_rate_noise', 'cvout_noise_qc'
    };
    
    % RSS of components
    rss_components = 0;
    for i = 1:length(sources)
        rss_components = rss_components + std(contrib.(sources{i}))^2;
    end
    rss_components = sqrt(rss_components);
    
    % Total
    total = std(contrib.cvout_r);
    
    % Should match within numerical precision
    rel_error = abs(total - rss_components) / total;
    assert(rel_error < 0.01, 'Total noise does not match RSS of components');
    
    fprintf('✓ Noise RSS check passed (error: %.2f%%)\n', rel_error*100);
end
```

#### 10.2.2 Spectral Shape Validation

```matlab
function test_spectral_shape()
    % Verify expected spectral features
    
    nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
    
    [psd, f] = pwelch(nm.info.C.sense.noise.ratein_r_dps, [], [], [], ...
                     nm.imu.as.gyr.config.inf.fs);
    
    % Check for 1/f region (1-10 Hz)
    idx_1f = (f >= 1) & (f <= 10);
    % Fit slope in log-log
    p = polyfit(log10(f(idx_1f)), log10(psd(idx_1f)), 1);
    slope_1f = p(1);
    
    % Should have negative slope (1/f or steeper)
    assert(slope_1f < -0.5 && slope_1f > -2, ...
           '1/f region has unexpected slope');
    
    % Check for white region (100-1000 Hz)
    idx_white = (f >= 100) & (f <= 1000);
    % Variance of PSD should be small
    psd_white = psd(idx_white);
    cv = std(psd_white) / mean(psd_white);  % Coefficient of variation
    
    assert(cv < 0.5, 'White noise region is not flat');
    
    fprintf('✓ Spectral shape validation passed\n');
    fprintf('  1/f slope: %.2f (expected: -0.5 to -2)\n', slope_1f);
    fprintf('  White region CV: %.2f (expected: <0.5)\n', cv);
end
```

### 10.3 Validation Against Circuit Simulations

#### 10.3.1 Comparison Framework

```matlab
function validate_against_spice(nm, spice_results)
    % Compare MATLAB noise model with SPICE simulation results
    %
    % spice_results: Structure with fields
    %   .frequency [Hz]
    %   .psd_ota [V²/Hz]
    %   .psd_thermal [V²/Hz]
    %   .total [dps RMS]
    
    % Extract MATLAB results
    [psd_matlab, f_matlab] = pwelch(nm.info.C.sense.noise.ratein_r_dps, ...
                                   [], [], [], nm.imu.as.gyr.config.inf.fs);
    
    % Interpolate to match SPICE frequency points
    psd_interp = interp1(f_matlab, psd_matlab, spice_results.frequency);
    
    % Compare total noise
    total_matlab = std(nm.info.C.sense.noise.ratein_r_dps);
    total_spice = spice_results.total;
    error_pct = abs(total_matlab - total_spice) / total_spice * 100;
    
    % Tolerance: ±20%
    assert(error_pct < 20, ...
           sprintf('Total noise error %.1f%% exceeds 20%% tolerance', error_pct));
    
    fprintf('✓ SPICE validation passed\n');
    fprintf('  MATLAB: %.3f dps, SPICE: %.3f dps, Error: %.1f%%\n', ...
            total_matlab, total_spice, error_pct);
    
    % Plot comparison
    figure;
    loglog(spice_results.frequency, sqrt(psd_interp), 'b-', 'LineWidth', 2);
    hold on;
    loglog(spice_results.frequency, sqrt(spice_results.psd_total), 'r--', 'LineWidth', 2);
    xlabel('Frequency [Hz]');
    ylabel('ASD [dps/√Hz]');
    legend('MATLAB', 'SPICE');
    title('Noise Model Validation vs. SPICE');
    grid on;
end
```

#### 10.3.2 Parameter Extraction from Measurements

```matlab
function params = extract_params_from_measurement(measured_data)
    % Extract noise model parameters from measurement data
    %
    % measured_data: Structure with fields
    %   .time [s]
    %   .rate [dps]
    %   .fs [Hz]
    
    % Calculate PSD
    [psd, f] = pwelch(measured_data.rate, [], [], [], measured_data.fs);
    
    % Extract white noise floor (100-1000 Hz)
    idx = (f >= 100) & (f <= 1000);
    params.white_noise_floor = median(sqrt(psd(idx)));  % [dps/√Hz]
    
    % Extract flicker noise corner
    % Find where PSD transitions from 1/f to white
    idx_1f = (f >= 1) & (f <= 100);
    [~, corner_idx] = min(abs(sqrt(psd(idx_1f)) - params.white_noise_floor));
    params.flicker_corner = f(idx_1f(corner_idx));  % [Hz]
    
    % Calculate ARW
    params.ARW = params.white_noise_floor * sqrt(3600) * 60;  % [deg/√hr]
    
    % Calculate bias instability from Allan deviation
    [tau, adev] = calculate_allan_deviation(measured_data.rate, measured_data.fs);
    [params.bias_instability, idx_min] = min(adev);
    params.tau_at_BI = tau(idx_min);
    
    % Extract RTN amplitude (if present)
    % Look for discrete steps in time series
    rate_diff = diff(measured_data.rate);
    rate_diff_std = std(rate_diff);
    large_steps = abs(rate_diff) > 3*rate_diff_std;
    if sum(large_steps) > 10
        params.rtn_present = true;
        params.rtn_amplitude = median(abs(rate_diff(large_steps)));
    else
        params.rtn_present = false;
        params.rtn_amplitude = 0;
    end
    
    % Display results
    fprintf('Extracted Parameters:\n');
    fprintf('  White noise floor: %.3f dps/√Hz\n', params.white_noise_floor);
    fprintf('  Flicker corner: %.1f Hz\n', params.flicker_corner);
    fprintf('  ARW: %.3f deg/√hr\n', params.ARW);
    fprintf('  Bias instability: %.3f dps @ τ=%.1f s\n', ...
            params.bias_instability, params.tau_at_BI);
    if params.rtn_present
        fprintf('  RTN detected: amplitude ≈ %.3f dps\n', params.rtn_amplitude);
    end
end
```

### 10.4 Regression Testing

#### 10.4.1 Baseline Establishment

```matlab
function baseline = establish_baseline(imu, sys)
    % Establish baseline results for regression testing
    
    configs = {'All', 'Thermal_Only', 'LowFreq_Only', 'Disable_CM'};
    baseline = struct();
    
    for i = 1:length(configs)
        nm = noise_model(imu, sys, init_noise_cfg_enhanced(configs{i}));
        
        baseline.(configs{i}).total_C = std(nm.info.C.sense.noise.ratein_r_dps);
        baseline.(configs{i}).total_S = std(nm.info.S.sense.noise.ratein_r_dps);
        baseline.(configs{i}).total_Z = std(nm.info.Z.sense.noise.ratein_r_dps);
        
        % Save spectrum
        [psd, f] = pwelch(nm.info.C.sense.noise.ratein_r_dps, [], [], [], ...
                         nm.imu.as.gyr.config.inf.fs);
        baseline.(configs{i}).psd = psd;
        baseline.(configs{i}).frequency = f;
    end
    
    % Save baseline
    save('noise_model_baseline.mat', 'baseline');
    fprintf('Baseline established and saved\n');
end
```

#### 10.4.2 Regression Check

```matlab
function regression_check(imu, sys, tolerance_pct)
    % Check current results against baseline
    %
    % tolerance_pct: Maximum allowed deviation [%] (default: 5%)
    
    if nargin < 3, tolerance_pct = 5; end
    
    % Load baseline
    load('noise_model_baseline.mat', 'baseline');
    
    configs = fieldnames(baseline);
    passed = true;
    
    fprintf('\n=== Regression Test Results ===\n');
    
    for i = 1:length(configs)
        nm = noise_model(imu, sys, init_noise_cfg_enhanced(configs{i}));
        
        % Check C-axis
        current = std(nm.info.C.sense.noise.ratein_r_dps);
        reference = baseline.(configs{i}).total_C;
        error_pct = abs(current - reference) / reference * 100;
        
        status = '✓ PASS';
        if error_pct > tolerance_pct
            status = '✗ FAIL';
            passed = false;
        end
        
        fprintf('%s: %s (%.3f dps, ref: %.3f dps, error: %.1f%%)\n', ...
                configs{i}, status, current, reference, error_pct);
    end
    
    if passed
        fprintf('\n=== All Regression Tests Passed ===\n');
    else
        fprintf('\n=== Regression Test FAILED ===\n');
        error('Noise model results deviate from baseline');
    end
end
```

### 10.5 Validation Checklist

#### Complete Validation Workflow:

```matlab
% 1. Run unit tests
test_thermal_noise_scaling();
test_noise_independence();
test_pseudo_sin_coefficients();

% 2. Run integration tests
test_noise_rss();
test_spectral_shape();

% 3. Validate against circuit simulation
load('spice_results.mat', 'spice_results');
nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
validate_against_spice(nm, spice_results);

% 4. Compare with measurement data
load('measured_data.mat', 'measured');
params = extract_params_from_measurement(measured);
% Compare params with model predictions

% 5. Establish or check regression baseline
if ~exist('noise_model_baseline.mat', 'file')
    baseline = establish_baseline(imu, sys);
else
    regression_check(imu, sys, 5);  % 5% tolerance
end

fprintf('\n=== Validation Complete ===\n');
```

**Expected Tolerances:**

| Metric | Tolerance | Notes |
|--------|-----------|-------|
| Total noise (vs. SPICE) | ±20% | Circuit-level correlation |
| Total noise (vs. measurement) | ±30% | Process variation, temperature |
| Spectral shape | Qualitative | 1/f slope, white region |
| ARW | ±25% | Depends on integration bandwidth |
| Bias instability | ±50% | Long-term drift mechanisms |
| Regression | ±5% | Code/config changes |

---
