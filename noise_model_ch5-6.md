
## 5. Sigma-Delta Modulator

### 5.1 MASH 1-1 Architecture Overview

The rate output uses a Multi-stAge noise SHaping (MASH) architecture with two first-order stages.

**Architecture:**

```
           ┌─────────────────────────────────────────┐
  Input    │  1st Order SD                           │  v1
  ────────>│  (3-level)    ────────────────────────> │────────┐
           │               Quantization              │        │
           │               Error y1                  │        │
           └────────────────┬────────────────────────┘        │
                            │                                 │
                            v                                 │
           ┌────────────────────────────────────────┐        │
  Error y1 │  2nd Order SD                          │  v2    │
  ────────>│  (2-level)    ────────────────────────> │────────┤
           │                                         │        │
           └─────────────────────────────────────────┘        │
                                                               │
           ┌─────────────────────────────────────────┐        │
           │  Digital Combiner                       │ Output │
           │  H(z) = 1 - z^(-1)                      │<───────┘
           │  Cancels 1st stage Q-noise              │
           └─────────────────────────────────────────┘
```

**Advantages:**
- Second-order noise shaping: NTF = $(1 - z^{-1})^2$
- Inherently stable (two 1st-order stages)
- Digital error cancellation
- Good SNR in signal bandwidth

### 5.2 First-Stage Sigma-Delta Modulator

#### 5.2.1 Input Signal Formation

The input to SD1 includes rate signal and all sense chain noise:

```matlab
% SD1 input formation
% File: sim_rate_sd_noise.m, Line 10

sd.C.u = obj.imu.as.gyr.config.inf.sense.cv.C.cv_gain * ...
         (obj.sys.rate_x.Value * cos(2*pi*ctrl.rate_freq*obj.ctrl.t) + ...
          obj.info.C.sense.noise.ratein_r_dps);
```

#### 5.2.2 Phase Noise Modulation of Feedback

**Critical Feature:** The SD1 feedback DAC switches at the drive frequency. Drive phase noise modulates the feedback timing, creating a multiplicative noise term.

**Modulation Factor:**

$$\gamma = \frac{2}{N_{levels,SD1}}$$

For 3-level SD: $\gamma = 2/3$

**Feedback Modulation:**

$$FB_{mod}(t) = 1 + \gamma \times \left[\phi_{drive,noise}(t) + \frac{I_{DAC,noise}(t)}{I_{DAC,nom}}\right]$$

**Implementation:**

```matlab
% Phase noise modulation of SD feedback
% File: sim_sense_noise.m, Line 214

sns_noise.pll.sd_mod = ones(1,obj.ctrl.N) + ...
    (obj.ctrl.noise_sel('Qnoise_dem') * ...
     (2/sense.rate.common.sd1_nlev) * ...
     (drive.noise.pll.noise.phase_noise_drive + ...
      (sqrt(2)*idac_bg_noise/ref.idac.iout)));
```

**Physical Mechanism:**
1. Phase jitter shifts DAC switching time
2. Creates charge transfer error
3. Modulates feedback level
4. **Folds high-frequency quantization noise into signal band**

This is a unique noise coupling mechanism in sigma-delta modulators with high-frequency feedback clocking.

#### 5.2.3 SD1 Simulation

**Modified DSM Simulator:**

```matlab
% Simulate 1st SD with phase noise modulation
% File: sim_rate_sd_noise.m, Line 22

[sd.C.v1, xn, xmax, sd.C.y1] = simulateDSM_mod(...
    sd.C.u, ...                          % Input signal + noise
    sd.ABCD_model_SD1, ...              % State-space model
    obj.imu.as.gyr.config.inf.sense.rate.common.sd1_nlev, ... % Levels
    [], ...                             % No input-referred dithering
    obj.info.sense.noise.pll.sd_mod);   % Phase noise modulation

% Normalize output: +2/0/-2 → +1/0/-1
sd.C.v1 = sd.C.v1 / 2;
```

**Outputs:**
- `v1`: Digital bitstream output [±1, 0]
- `y1`: Quantization error (input to SD2)

#### 5.2.4 Quantization Noise Extraction

The quantization error is extracted for the second stage:

```matlab
% Calculate SD1 quantization noise
% File: sim_rate_sd_noise.m, Line 31

% Use differencing filter to extract error
sd.C.u2 = lsim(sd.Hsd_model_diff, [sd.C.v1 ; sd.C.y1], obj.ctrl.t)';
```

### 5.3 Second-Stage Sigma-Delta Modulator

#### 5.3.1 Purpose

The second stage:
1. Shapes the quantization noise from SD1
2. Provides second-order overall noise shaping
3. Enables digital cancellation of SD1 quantization noise

#### 5.3.2 SD2 Simulation

```matlab
% Simulate 2nd SD
% File: sim_rate_sd_noise.m, Line 51

[sd.C.v2, xn2, xmax2, y2] = simulateDSM(...
    sd.C.u2, ...              % SD1 quantization error
    sd.ABCD_model_SD2, ...    % State-space model
    2, ...                    % 2-level quantizer
    []);                      % No phase noise modulation (internal clock)
```

**Note:** SD2 uses an internal clock (not drive frequency), so no phase noise modulation.

### 5.4 Digital Combiner

#### 5.4.1 Error Cancellation Filter

The digital combiner cancels the first-stage quantization noise:

$$H_{dig}(z) = 1 - z^{-1}$$

This differentiator:
- Cancels SD1 quantization noise
- Passes SD2 output with first-order shaping
- Results in second-order shaping overall

**Implementation:**

```matlab
% Digital combination
% File: sim_rate_sd_noise.m, Line 72

sd.C.v = lsim(sd.Hsd_model_dig, [sd.C.v1 ; sd.C.v2], obj.ctrl.t)';
```

#### 5.4.2 Overall Transfer Functions

**Signal Transfer Function (STF):**

$$H_{STF}(z) = z^{-1}$$

(One sample delay)

**Noise Transfer Function (NTF):**

$$H_{NTF}(z) = (1 - z^{-1})^2$$

(Second-order high-pass)

**Calculation:**

```matlab
% Calculate transfer functions
[sd.C.ntf, sd.C.stf] = calculateTF(sd.ABCD_model_SD1, sd.C.gain_quantizer1);
```

### 5.5 Quantization Noise Shaping Visualization

**Typical NTF Shape:**

```
|NTF(f)|²
  ^
  |        
  |                         ___---
  |                   __----
  |             __----
  |       __----
  | __----    Second-order slope: +40 dB/decade
  |/
  +---------------------------------> f
  0              f_sig              fs/2
             (0-200 Hz)         (50 kHz)
```

**In Signal Band (0-200 Hz):**
- Quantization noise heavily suppressed
- Thermal/flicker noise dominant

**At High Frequencies (10-50 kHz):**
- Quantization noise rises rapidly
- Decimation filter removes this
- Only shaped noise in signal band remains

### 5.6 Practical Considerations

#### 5.6.1 OSR (Oversampling Ratio)

$$OSR = \frac{f_s}{2 \times BW_{signal}} = \frac{100 \text{ kHz}}{2 \times 200 \text{ Hz}} = 250$$

High OSR provides:
- More quantization noise suppression
- Better SNR
- Allows simpler anti-aliasing filter

#### 5.6.2 Stability

MASH architecture is inherently stable because:
- Each stage is first-order (unconditionally stable)
- No feedback between stages
- Digital cancellation in post-processing

#### 5.6.3 Decimation Filtering

After the sigma-delta modulator, a digital decimation filter:
1. Low-pass filters to remove shaped quantization noise
2. Decimates from 100 kHz to final output rate (typically 100-1000 Hz)
3. Provides additional anti-aliasing

**Not explicitly modeled in this noise analysis** (occurs in digital domain after all analog noise is shaped).

---

## 6. Pseudo-Sinusoidal Demodulation

### 6.1 Motivation

Pure sinusoidal multiplication for demodulation is difficult to implement in switched-capacitor circuits. A **pseudo-sinusoidal** approximation uses:
- Switched resistors (simple to implement)
- Step-wise approximation to sine wave
- Digital-like control signals

### 6.2 Waveform Definition

#### 6.2.1 Period and Levels

**Period:** 192 samples (at 100 kHz sampling → ~520 Hz pseudo-frequency)

**Amplitude Levels:** Three levels plus zero
- High: 1.0
- Medium: 1/2.4 ≈ 0.417
- Low: -1.0
- Medium-low: -1/2.4 ≈ -0.417
- Zero: 0

#### 6.2.2 Phase Assignments

```matlab
% Pseudo-sinusoidal waveform generator
% File: noise_model.m (methods section)

function w = pseudo_sin(t, N)
    % N: Period in samples (192 typical)
    % t: Time indices (0 to length-1)
    
    frac = t ./ N;  % Normalized phase [0, 1)
    w = zeros(size(t));
    
    % Positive half-cycle
    w(frac >= 4/192 & frac < 24/192) = 1/2.4;    % Rising edge
    w(frac >= 24/192 & frac < 72/192) = 1.0;     % Peak
    w(frac >= 72/192 & frac < 94/192) = 1/2.4;   % Falling edge
    
    % Zero crossing
    w(frac >= 94/192 & frac < 98/192) = 0;
    
    % Negative half-cycle
    w(frac >= 98/192 & frac < 120/192) = -1/2.4;  % Falling edge
    w(frac >= 120/192 & frac < 168/192) = -1.0;   % Trough
    w(frac >= 168/192 & frac < 190/192) = -1/2.4; % Rising edge
    
    % Zero crossing (0-4 and 190-192)
    % Already initialized to 0
end
```

#### 6.2.3 Visual Representation

```
 1.0  ────────────────────────
 0.4  ────                    ────
 0.0  ────                        ────────
-0.4                                     ────        ────
-1.0                                         ────────────
     |  |           |        |  |  |      |           |  |
     0  4          24       72 94 98    120         168 190
                                Phase [samples]
```

Compare to true sine wave:
- Approximates well near zero crossings
- Some harmonic distortion at peaks
- **Total Harmonic Distortion (THD) ≈ 15-20%**

### 6.3 Signal and Noise Coefficients

#### 6.3.1 Signal Coefficient $K_{in}$

The signal coefficient relates integrated charge to input signal:

$$Q_{integrated} = K_{in} \times V_{input}$$

**Calculation:**

$$K_{in} = \int_0^T V_{signal}(t) \times w_{pseudo}(t) \, dt$$

where $w_{pseudo}(t)$ is the pseudo-sin waveform.

**From Equations (Section 2.3.2):**

For pseudo-sin with phase integrator:
- Depends on switched resistor values
- Includes effects of timing phases ($\phi_{hi}$, $\phi_{lo}$, $\phi_{grd}$)
- Extracted through circuit simulation or measurement

**Typical Value:** $K_{in} \approx 0.7-0.9$ (normalized)

#### 6.3.2 Noise Coefficient $K_{in,noise}$

The noise coefficient differs from signal coefficient because noise is:
- Non-coherent with demodulation
- Averaged differently by switching

**Calculation:**

$$K_{in,noise} = \sqrt{\int_0^T w_{pseudo}^2(t) \, dt}$$

This is the RMS value of the pseudo-sin waveform.

**Theoretical Calculation:**

For the waveform defined above:
```matlab
% Calculate noise coefficient
w = pseudo_sin(0:191, 192);
K_in_noise_theory = sqrt(mean(w.^2));
% Result: K_in_noise ≈ 0.7613
```

**Practical Value:** $K_{in,noise} \approx 0.761$ (from FFT analysis)

#### 6.3.3 Ratio and Physical Meaning

**Coefficient Ratio:**

$$\frac{K_{in,noise}}{K_{in}} \approx \frac{0.761}{0.85} \approx 0.895$$

**Physical Interpretation:**
- Signal experiences full correlation over period → higher coefficient
- Noise averages with RMS weighting → lower coefficient
- Ratio represents "noise bandwidth penalty" of pseudo-sin vs. true sine

### 6.4 Application in Noise Model

#### 6.4.1 Signal Scaling

When scaling signal from voltage to rate:

```matlab
% Signal uses K_in
signal_rate = (signal_voltage / cv_gain) * (FSdps / Vref) * ...
              sense.rate.coeff.Kin;
```

#### 6.4.2 Noise Scaling

When transferring noise through demodulator:

```matlab
% Noise uses both coefficients
noise_at_integrator_input = noise_at_cv_output * ...
                           (sense.rate.coeff.Kin_noise / sense.rate.coeff.Kin);
```

**Example from Code:**

```matlab
% Total noise at Rate Int Input
% File: sim_sense_noise.m, Line 185

C_sense.ratein_r = C_sense.cvout_r * ...
                   sense.rate.coeff.Kin_noise / sense.rate.coeff.Kin;
```

#### 6.4.3 Reference Noise

Reference noise (bandgap, IDAC) goes through different path:

```matlab
% Reference noise coefficient ratio
(sense.rate.coeff.Kref_noise / sense.rate.coeff.Kin_noise)
```

This ratio accounts for how reference and signal demodulate differently through the pseudo-sin switching.

### 6.5 Pseudo-Sin Modulation of Colored Noise

For the sense C/V OTA, colored noise is explicitly modulated:

```matlab
% Generate pseudo-sin waveform for full signal length
pseudo_sin_signal = repmat(obj.pseudo_sin(0:ceil(obj.ctrl.N/192)*192-1, 192), ...
                           ceil(obj.ctrl.N/192));

% Apply modulation
C_sense.cvout_cvs_ota = ...
    (sns_noise.cv.noise.cvs_noise * sns_noise.cv.noise.C.cvs_gain) .* ...
    pseudo_sin_signal(1:obj.ctrl.N);
```

**Why This Matters:**
- Translates baseband noise to drive frequency
- Creates sidebands around $f_d$
- Only properly demodulated components reach output
- Effectively band-limits the noise

### 6.6 Design Trade-offs

#### 6.6.1 Pseudo-Sin vs. True Sine

**Advantages of Pseudo-Sin:**
- Simple switched-capacitor implementation
- No analog multipliers needed
- Digital control signals (easy to generate)
- Good enough for moderate performance gyros

**Disadvantages:**
- Harmonic distortion (15-20% THD)
- Noise bandwidth penalty (~10%)
- Slightly reduced SNR vs. true sine

#### 6.6.2 Optimization

**Key Parameter:** Number of levels and phase assignments

Current implementation (3 levels) balances:
- Circuit complexity (fewer switches)
- THD performance
- Noise coefficient

**Alternative:** 5-level or 7-level pseudo-sin
- Better sine approximation
- Lower THD
- More complex switching
- Diminishing returns beyond 5 levels

---
