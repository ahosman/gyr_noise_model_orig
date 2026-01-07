## 7. Random Telegraph Noise (RTN)

### 7.1 Physical Origin

Random Telegraph Noise arises from charge trapping and detrapping in semiconductor oxide traps, particularly in MOSFET gate oxides.

**Physical Process:**

```
Trap State Transitions:
         
  Empty  ──k_capture──>  Filled
  State  <──k_emit────   State
         
  k_capture = f(trap energy, temperature, bias)
  k_emit    = f(trap energy, temperature)
```

**Characteristics:**
- **Discrete amplitude steps** (not continuous like thermal noise)
- **Random timing** (Poisson statistics)
- **Temperature dependent** (Arrhenius behavior)
- **Bias dependent** (trap occupancy)

### 7.2 Mathematical Model

#### 7.2.1 Two-State Markov Process

RTN from a single trap is modeled as a telegraph signal:

$$V_{RTN}(t) = \begin{cases}
+A & \text{trap filled (electron captured)} \\
0 & \text{trap empty (electron emitted)}
\end{cases}$$

**Time Constants:**
- $\tau_{high}$: Mean time in filled state [s]
- $\tau_{low}$: Mean time in empty state [s]

**Transition Probabilities:**

Per sample (time step $\Delta t$):

$$P_{high \to low} = 1 - e^{-\Delta t / \tau_{high}} \approx \frac{\Delta t}{\tau_{high}} \quad \text{for } \Delta t \ll \tau_{high}$$

$$P_{low \to high} = 1 - e^{-\Delta t / \tau_{low}} \approx \frac{\Delta t}{\tau_{low}} \quad \text{for } \Delta t \ll \tau_{low}$$

#### 7.2.2 Power Spectral Density

Single trap generates a Lorentzian spectrum:

$$S_{RTN}(f) = \frac{4 A^2 \tau_{avg}}{1 + (2\pi f \tau_{avg})^2}$$

where:

$$\tau_{avg} = \frac{\tau_{high} \times \tau_{low}}{\tau_{high} + \tau_{low}}$$

**Corner Frequency:**

$$f_c = \frac{1}{2\pi \tau_{avg}}$$

**Spectrum Shape:**
- White (flat) for $f \ll f_c$
- Roll-off at 20 dB/decade for $f \gg f_c$

#### 7.2.3 Multiple Traps

Real transistors have multiple traps with distributed time constants:

$$S_{RTN,total}(f) = \sum_{i=1}^{N_{traps}} \frac{4 A_i^2 \tau_{avg,i}}{1 + (2\pi f \tau_{avg,i})^2}$$

**Log-Distributed Traps:**

If $\tau_i$ spans several decades (e.g., 1 μs to 1 ms):
- Creates approximately **1/f spectrum** at low frequencies
- White plateau at high frequencies
- Looks similar to flicker noise but different origin

### 7.3 Implementation

#### 7.3.1 RTN Generation Function

```matlab
function rtn_signal = generate_rtn(N, fs, amplitude, tau_high, tau_low, num_traps, seed)
% GENERATE_RTN - Generate Random Telegraph Noise signal
%
% Inputs:
%   N         - Number of samples
%   fs        - Sampling frequency [Hz]
%   amplitude - RTN step amplitude [V] (scalar or vector)
%   tau_high  - Mean time in trapped state [s] (scalar or vector)
%   tau_low   - Mean time in untrapped state [s] (scalar or vector)
%   num_traps - Number of independent traps (default: 1)
%   seed      - Random seed for reproducibility (default: random)
%
% Output:
%   rtn_signal - RTN time series [V]

    if nargin < 6, num_traps = 1; end
    if nargin < 7, rng('shuffle'); else rng(seed); end
    
    % Expand scalars to vectors
    if isscalar(amplitude), amplitude = amplitude * ones(1, num_traps); end
    if isscalar(tau_high), tau_high = tau_high * ones(1, num_traps); end
    if isscalar(tau_low), tau_low = tau_low * ones(1, num_traps); end
    
    % Initialize
    rtn_signal = zeros(1, N);
    dt = 1/fs;
    
    % Generate RTN for each trap
    for trap = 1:num_traps
        A = amplitude(trap);
        tau_h = tau_high(trap);
        tau_l = tau_low(trap);
        
        % Transition probabilities
        p_high_to_low = min(1 - exp(-dt/tau_h), 0.5);
        p_low_to_high = min(1 - exp(-dt/tau_l), 0.5);
        
        % Initial state (equilibrium probability)
        state = rand() > (tau_l/(tau_h + tau_l));
        trap_signal = zeros(1, N);
        
        % Markov chain simulation
        for i = 1:N
            trap_signal(i) = A * double(state);
            
            % State transition
            if state
                if rand() < p_high_to_low, state = 0; end
            else
                if rand() < p_low_to_high, state = 1; end
            end
        end
        
        % Remove DC and add to total
        rtn_signal = rtn_signal + (trap_signal - mean(trap_signal));
    end
end
```

#### 7.3.2 Integration into Noise Model

**Where to Apply RTN:**

1. **Sense C/V OTA Input Transistors** (most critical)
2. **IDAC Current Source Transistors**
3. **VCO Bias Circuit** (affects phase noise)
4. **Phase Integrator OTA**

**Example Integration:**

```matlab
% Add RTN to Sense C/V OTA
% Enhancement to sim_sense_noise.m

if obj.ctrl.noise_sel('CVS_RTN')
    % Define trap parameters (from characterization)
    num_traps = 3;
    tau_high = [15e-6, 30e-6, 80e-6];  % Time constants [s]
    tau_low = [12e-6, 25e-6, 70e-6];
    amplitude = [1.5e-3, 2.0e-3, 1.2e-3];  % Voltage steps [V]
    
    % Generate RTN signal
    rtn_cvs = generate_rtn(obj.ctrl.N, obj.imu.as.gyr.config.inf.fs, ...
                          amplitude, tau_high, tau_low, num_traps, seed_cvs);
    
    % Add to C/V noise with appropriate gain
    C_sense.cvout_cvs_rtn = rtn_cvs * sns_noise.cv.noise.C.cvs_gain;
    
    % Include in total
    C_sense.cvout_cvs_ota = C_sense.cvout_cvs_ota + C_sense.cvout_cvs_rtn;
end
```

### 7.4 Temperature Dependence

RTN time constants follow Arrhenius behavior:

$$\tau(T) = \tau_0 \times \exp\left(\frac{E_a}{k_B T}\right)$$

where:
- $E_a$: Activation energy (trap-dependent, typically 0.3-1.0 eV)
- $k_B$: Boltzmann constant
- $T$: Absolute temperature [K]

**Temperature Scaling:**

```matlab
function tau_scaled = scale_rtn_tau(tau_ref, T_ref, T_current, Ea)
% Scale RTN time constant with temperature
%
% Ea: Activation energy [eV] (typical: 0.5 eV)
    kB_eV = 8.617e-5;  % Boltzmann constant [eV/K]
    tau_scaled = tau_ref * exp(Ea/kB_eV * (1/T_current - 1/T_ref));
end

% Example: Scale from 25°C to 85°C
tau_25C = 20e-6;  % 20 μs at 25°C
Ea = 0.5;  % 0.5 eV activation energy
tau_85C = scale_rtn_tau(tau_25C, 273.15+25, 273.15+85, Ea);
% Result: tau_85C ≈ 5 μs (faster at higher temperature)
```

### 7.5 Characterization Methods

#### 7.5.1 Time-Domain Analysis

**Identify Traps:**
1. Record long time series (>1 second)
2. Look for discrete amplitude steps
3. Measure step amplitudes ($A_i$)
4. Measure dwell times (histogram → $\tau_i$)

#### 7.5.2 Frequency-Domain Analysis

**Extract Parameters from PSD:**

1. Measure low-frequency PSD
2. Fit Lorentzian functions:
   $$S(f) = \sum_i \frac{4A_i^2 \tau_{avg,i}}{1 + (2\pi f \tau_{avg,i})^2}$$
3. Extract $A_i$ and $\tau_{avg,i}$ from fits
4. Validate in time domain

#### 7.5.3 Statistical Analysis

**From Multiple Devices:**

- Amplitude distribution (log-normal typical)
- Time constant distribution (very wide, log-uniform)
- Trap density (number per device)
- Temperature dependence ($E_a$ extraction)

### 7.6 Design Guidelines

#### 7.6.1 When RTN Matters

RTN is significant when:
- **Transistor area < 1 μm²** (W×L)
- **Operating frequency < 1/τ_avg** (trap activity visible)
- **PMOS input pairs** (p+ poly has more traps)
- **Temperature > 85°C** (faster trap dynamics)
- **Input-referred** (OTA, comparator inputs)

#### 7.6.2 Mitigation Strategies

1. **Larger Input Transistors**
   - W×L > 10 μm² minimizes effect
   - Dilutes single-trap impact
   - Trade-off: larger area, more capacitance

2. **Parallel Transistors**
   - Use M parallel devices (M=4-8)
   - Averaging reduces peak amplitude by √M
   - Uncorrelated trap events

3. **NMOS Input Pairs**
   - n+ poly gate has fewer interface traps
   - Not always possible (depends on circuit)

4. **Chopper Stabilization**
   - Modulate signal above RTN corner frequency
   - Demodulate after amplification
   - Requires $f_{chop} > 1/(2\pi\tau_{avg})$

5. **Process Improvements**
   - High-quality gate oxide growth
   - Post-oxidation annealing
   - Thicker oxide (if voltage allows)

### 7.7 Typical Parameters

**Consumer MEMS Gyroscope:**
- OTA input transistor: W/L = 20/1 μm → Area ≈ 20 μm²
- Number of traps: 2-5 dominant traps
- Amplitude: 0.5-3 mV (input-referred)
- Time constants: 1 μs - 1 ms
- Corner frequencies: 160 Hz - 160 kHz

**High-Performance Gyroscope:**
- Larger transistors (W×L > 50 μm²)
- Parallel input pairs (M=4)
- RTN contribution < 5% of total noise

---

## 8. Noise Budgeting and Analysis

### 8.1 Total Noise Calculation

#### 8.1.1 Sense Channel Noise at C/V Output

**Combining Independent Sources:**

All uncorrelated noise sources combine in RSS (Root Sum Square):

```matlab
% Generate final noise vector at C/V output [V]
% File: sim_sense_noise.m, Line 180

C_sense.cvout_r = ...
    C_sense.cvout_bg_noise +            % Bandgap + IDAC
    C_sense.cvout_r_cpcn +              % MEMS series R
    C_sense.cvout_cvs_ota +             % Sense C/V OTA
    C_sense.cvout_rsig +                % Rate int signal R
    C_sense.cvout_rate_ota +            % Rate int OTA
    C_sense.cvout_cm_c0_noise +         % CM over C0 offset
    C_sense.cvout_rp +                  % ASIC parallel R
    C_sense.cvout_cm_qdem_noise_r +     % CM via quad demod
    C_sense.cvout_cm_cpar_noise +       % CM via parasitic
    C_sense.cvout_cm_rate_noise +       % CM via rate signal
    C_sense.cvout_noise_qc;             % Quad comp coupling
```

**Note:** These are added as time-domain signals, each independently generated with `randn()` or colored noise generators using different seeds.

#### 8.1.2 Conversion to Rate Units

**At Rate Integrator Input:**

Accounting for pseudo-sin demodulation coefficients:

```matlab
% Noise at Rate Int Input [V]
C_sense.ratein_r = C_sense.cvout_r * ...
                   sense.rate.coeff.Kin_noise / sense.rate.coeff.Kin;

% Convert to rate units [dps]
C_sense.ratein_r_dps = (C_sense.ratein_r / sense.cv.C.cv_gain) + ...
                       C_sense.brw_noise_dps + ...
                       C_sense.quad_noise_dps;
```

**Components:**
1. **Electronic noise** (from C/V output) scaled by demod coefficients
2. **Brownian noise** (already in dps)
3. **Quadrature × phase noise** (already in dps)

#### 8.1.3 After Sigma-Delta Modulation

The SD modulator adds quantization noise (shaped):

```matlab
% SD output spectrum
spectrum_C_rate.sd.rate_out = plotFunction(...
    sd.C.v .* obj.imu.as.gyr.config.inf.sense.rate.C.FS, ...
    'fs', obj.imu.as.gyr.config.inf.fs, ...
    'f0', 0, 'fb', obj.ctrl.fb, 'fsig', ctrl.rate_freq, ...
    'stats', true, 'unit', 'dps', ...
    'plot_fft', obj.ctrl.plot_spectra_s, 'new_fig', obj.ctrl.plot_fft);
```

### 8.2 Noise Breakdown Reporting

#### 8.2.1 Automated Analysis

Use the `generate_noise_breakdown()` function:

```matlab
% Generate comprehensive noise report
opts.signal_band = 200;  % Hz
opts.plot = true;
opts.verbose = true;

report = generate_noise_breakdown(obj, opts);

% Access results
fprintf('Total Noise:\n');
fprintf('  C-axis: %.3f dps RMS\n', report.total.C);
fprintf('  S-axis: %.3f dps RMS\n', report.total.S);
fprintf('  Z-axis: %.3f dps RMS\n', report.total.Z);

fprintf('\nDominant Sources:\n');
fprintf('  C-axis: %s (%.1f%%)\n', ...
        report.dominant.C.source, report.dominant.C.percentage);
```

#### 8.2.2 Example Output

```
═══════════════════════════════════════════════════════════════════════════
 C (Roll) AXIS
═══════════════════════════════════════════════════════════════════════════

TOTAL NOISE: 0.0845 dps RMS (84.5 mdps RMS)

NOISE BREAKDOWN BY SOURCE:
───────────────────────────────────────────────────────────────────────────
Source                          RMS [mdps]   Contribution
───────────────────────────────────────────────────────────────────────────
Sense Cv Ota                        42.3      50.1% ████████████████████████
Brownian                            28.7      33.9% ████████████████
Mems Resistor                        8.2       9.7% ████
Rate Int Ota                         3.1       3.7% █
Cm C0 Offset                         1.8       2.1% █
Quad Comp                            0.9       1.1%
...
───────────────────────────────────────────────────────────────────────────
TOTAL                               84.5     100.0%
═══════════════════════════════════════════════════════════════════════════

DOMINANT SOURCE: Sense Cv Ota (50.1% of total)

FREQUENCY BAND ANALYSIS:
───────────────────────────────────────────────────────────────────────────
  DC-1Hz              :   12.345 mdps RMS
  1-10Hz              :   18.234 mdps RMS
  10-100Hz            :   38.456 mdps RMS
  100-200Hz           :   15.510 mdps RMS

KEY METRICS:
───────────────────────────────────────────────────────────────────────────
  Angle Random Walk (ARW):        0.085 deg/√hr
  Bias Instability (approx):     18.234 mdps
  Signal-to-Noise Ratio:         61.43 dB
  Effective Bits (ENOB):          9.9 bits
```

### 8.3 Frequency Domain Analysis

#### 8.3.1 Power Spectral Density

**Welch's Method:**

```matlab
% Calculate PSD
[psd, f] = pwelch(C_sense.ratein_r_dps, ...
                 hann(min(obj.ctrl.N/4, 8192)), ... % Window
                 [], ...                             % 50% overlap
                 [], ...                             % NFFT = window length
                 obj.imu.as.gyr.config.inf.fs);     % Sampling frequency

% Plot
figure;
loglog(f, sqrt(psd));  % ASD [dps/√Hz]
xlabel('Frequency [Hz]');
ylabel('Amplitude Spectral Density [dps/√Hz]');
grid on;
```

**Expected Features:**
- 1/f region: 1-100 Hz (flicker noise, RTN)
- White region: 100 Hz - 10 kHz (thermal noise dominant)
- Roll-off: > 10 kHz (signal bandwidth limit)

#### 8.3.2 Integrated Noise in Bands

**Calculate RMS in Specific Bands:**

```matlab
function rms_band = integrate_psd_band(psd, f, f_low, f_high)
    idx = (f >= f_low) & (f <= f_high);
    power_band = sum(psd(idx)) * mean(diff(f));  % Integrate PSD
    rms_band = sqrt(power_band);
end

% Example
rms_1_10Hz = integrate_psd_band(psd, f, 1, 10);
rms_10_100Hz = integrate_psd_band(psd, f, 10, 100);
rms_100_200Hz = integrate_psd_band(psd, f, 100, 200);
```

### 8.4 Allan Deviation

#### 8.4.1 Calculation

Allan deviation characterizes long-term stability:

```matlab
function [tau, adev] = calculate_allan_deviation(rate_signal, fs)
    % rate_signal: Rate time series [dps]
    % fs: Sampling frequency [Hz]
    
    N = length(rate_signal);
    tau = logspace(log10(1/fs), log10(N/fs/10), 50);
    adev = zeros(size(tau));
    
    for i = 1:length(tau)
        m = round(tau(i) * fs);  % Samples per average
        if m < 1, m = 1; end
        if m > N/2, break; end
        
        % Reshape and calculate
        n_groups = floor(N/m);
        reshaped = reshape(rate_signal(1:n_groups*m), m, n_groups);
        averages = mean(reshaped, 1);
        
        % Allan deviation formula
        adev(i) = std(diff(averages)) / sqrt(2);
    end
    
    % Truncate to valid points
    adev = adev(1:i);
    tau = tau(1:i);
end
```

#### 8.4.2 Interpretation

**Allan Deviation Plot Features:**

```
σ(τ)
  ^
  |  \
  |   \              Minimum
  |    \___       (Bias Instability)
  |        \___/
  |            \
  |             \
  +-------------------> τ
     ARW      BI      Rate Random Walk
   Region   Region      Region
```

**Key Metrics:**

1. **Angle Random Walk (ARW):** Slope = -1/2 at short τ
   $$ARW = \sigma(\tau=1) \times 60 \quad [\text{deg}/\sqrt{\text{hr}}]$$

2. **Bias Instability (BI):** Minimum of Allan deviation
   - Typically at τ = 10-100 seconds
   - Units: dps or deg/hr

3. **Rate Random Walk (RRW):** Slope = +1/2 at long τ
   - Due to drift mechanisms
   - Less critical for MEMS gyroscopes

### 8.5 Key Performance Metrics

#### 8.5.1 Angle Random Walk (ARW)

**Definition:** White noise floor expressed per square root of integration time.

**Calculation:**

```matlab
% From flat region of ASD (10-100 Hz typical)
idx_flat = (f >= 10) & (f <= 100);
white_noise_floor = median(sqrt(psd(idx_flat)));  % [dps/√Hz]

% Convert to deg/√hr
ARW_deg_sqrt_hr = white_noise_floor * sqrt(3600) * 60;
```

**Typical Values:**
- Consumer: 0.1-1.0 deg/√hr
- Automotive: 0.01-0.1 deg/√hr
- Tactical: 0.001-0.01 deg/√hr
- Navigation: < 0.001 deg/√hr

#### 8.5.2 Bias Instability

**From Allan Deviation:**

```matlab
[tau_allan, adev_allan] = calculate_allan_deviation(rate_signal, fs);
[BI_dps, min_idx] = min(adev_allan);
tau_at_min = tau_allan(min_idx);

fprintf('Bias Instability: %.3f dps at τ = %.1f s\n', BI_dps, tau_at_min);
```

**Typical Values:**
- Consumer: 10-100 dps
- Automotive: 1-10 dps
- Tactical: 0.1-1 dps
- Navigation: < 0.1 dps

#### 8.5.3 Signal-to-Noise Ratio (SNR)

**For Full-Scale Range:**

```matlab
FS_dps = 100;  % Full scale: ±100 dps typical
noise_rms = std(rate_signal);  % RMS noise [dps]

SNR_dB = 20*log10(FS_dps / noise_rms);
```

#### 8.5.4 Effective Number of Bits (ENOB)

**From SNR:**

$$ENOB = \frac{SNR_{dB} - 1.76}{6.02}$$

```matlab
ENOB = (SNR_dB - 1.76) / 6.02;
```

**Typical Values:**
- 8-10 bits: Consumer gyroscope
- 10-12 bits: Automotive/industrial
- 12-14 bits: Tactical grade
- >14 bits: Navigation grade

### 8.6 Sensitivity Analysis

#### 8.6.1 Parameter Sweep

**Identify Critical Parameters:**

```matlab
% Example: Sweep sense C/V OTA noise
vn_opa_range = linspace(5e-9, 15e-9, 20);  % 5-15 nV/√Hz
total_noise = zeros(size(vn_opa_range));

for i = 1:length(vn_opa_range)
    % Modify parameter
    imu.as.gyr.config.inf.sense.cv.noise.Vn_opa = vn_opa_range(i);
    
    % Run simulation
    nm = noise_model(imu, sys, 'All');
    
    % Extract total noise
    total_noise(i) = std(nm.info.C.sense.noise.ratein_r_dps);
end

% Plot sensitivity
figure;
plot(vn_opa_range*1e9, total_noise*1000, 'LineWidth', 2);
xlabel('Sense C/V OTA Noise [nV/√Hz]');
ylabel('Total Noise [mdps RMS]');
grid on;
```

#### 8.6.2 Correlation Analysis

**Between Parameters and Total Noise:**

```matlab
% Monte Carlo with parameter variations
N_runs = 100;
parameters = zeros(N_runs, 5);  % 5 key parameters
total_noise_mc = zeros(N_runs, 1);

for i = 1:N_runs
    % Vary parameters (example: ±20% uniform)
    imu.as.gyr.config.inf.sense.cv.noise.Vn_opa = ...
        Vn_opa_nom * (0.8 + 0.4*rand());
    % ... vary other parameters
    
    % Run simulation
    nm = noise_model(imu, sys, 'All');
    total_noise_mc(i) = std(nm.info.C.sense.noise.ratein_r_dps);
    
    % Store parameters
    parameters(i,1) = imu.as.gyr.config.inf.sense.cv.noise.Vn_opa;
    % ... store other parameters
end

% Calculate correlation
R = corrcoef([parameters, total_noise_mc]);
correlations = R(end, 1:end-1);

% Display
bar(correlations);
set(gca, 'XTickLabel', {'OTA Noise', 'Rs MEMS', 'CM Noise', ...
                        'Phase Noise', 'C0 Offset'});
ylabel('Correlation with Total Noise');
title('Parameter Sensitivity Analysis');
```

### 8.7 Design Optimization

#### 8.7.1 Noise Budget Allocation

**Top-Down Approach:**

1. **System Requirement:** Total noise < 0.1 dps RMS (example)

2. **Allocate to Subsystems:**
   - MEMS (Brownian): 0.03 dps → 9% of budget
   - Sense C/V: 0.05 dps → 25% of budget
   - Rate Integrator: 0.02 dps → 4% of budget
   - Phase Noise: 0.01 dps → 1% of budget
   - CM Coupling: 0.02 dps → 4% of budget
   - Other: 0.01 dps → 1% of budget
   - **Total RSS:** $\sqrt{0.03^2 + 0.05^2 + ... } = 0.074$ dps
   - **Margin:** 35% (accounts for correlation, unknowns)

3. **Component-Level Specs:**
   - From Sense C/V budget → OTA noise < 8 nV/√Hz
   - From MEMS budget → Q > 10000, m > 0.5 mg
   - Etc.

#### 8.7.2 Optimization Strategy

**Identify Dominant Sources:**

```matlab
report = generate_noise_breakdown(nm);

% Sort by contribution
[~, idx] = sort(structfun(@(x) x, rmfield(report.breakdown.C, 'total')), 'descend');
source_names = fieldnames(rmfield(report.breakdown.C, 'total'));

fprintf('Optimization Priority:\n');
for i = 1:5  % Top 5 sources
    fprintf('%d. %s (%.1f%% contribution)\n', i, ...
            source_names{idx(i)}, report.percentage.C.(source_names{idx(i)}));
end
```

**Focus Effort on Top 2-3 Contributors:**
- 80/20 rule: Top 20% of sources contribute 80% of noise
- Optimize high-impact parameters first
- Diminishing returns beyond top 3 sources

---
