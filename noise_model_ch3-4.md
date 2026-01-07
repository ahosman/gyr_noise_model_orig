
## 3. Sense Channel Noise

### 3.1 Sense Channel Architecture

The sense channel detects Coriolis-induced capacitance changes and converts them to a digital rate output.

```
MEMS Sense    Sense C/V    Rate          Sigma-Delta    Digital
Electrodes → Converter → Integrator →   Modulator   →  Output
 (CP/CN)       (CVS)       (Rint)       (MASH 1-1)      (Rate)
                ↑            ↑              ↑
             CM Noise    Quad Comp     Phase Noise
             Coupling    Coupling      SD Modulation
```

**Key Features:**
- Fully differential architecture (CP/CN electrodes)
- Pseudo-sinusoidal demodulation at drive frequency
- Rate feedback via sigma-delta DAC
- Three independent axes (C, S, Z)

### 3.2 Reference Voltage Noise Sources

#### 3.2.1 Bandgap Reference

The bandgap voltage reference generates the system voltage reference (~1.2V typical). It exhibits both thermal and flicker noise.

**Colored Noise Generation:**

```matlab
% Bandgap reference noise
% File: sim_sense_noise.m, Line 26

[~, ref.bg.noise.bg_noise] = colouredNoise(...
    (obj.ctrl.N-1)*1/obj.imu.as.gyr.config.inf.fs, ...  % Duration
    obj.imu.as.gyr.config.inf.fs, ...                   % Sampling freq
    ref.bg.noise.freq, ...                              % Frequency array
    obj.ctrl.noise_sel('Vref') * ref.bg.noise.vnd, ... % PSD array
    43526, 0, 0, 1, false);                             % Seed, options
```

**Typical Spectrum:**
- White noise floor: 50-100 nV/√Hz
- Flicker corner: 100 Hz - 1 kHz
- 1/f region: increases at 10 dB/decade below corner

#### 3.2.2 IDAC Current Source

The IDAC generates the reference current for the rate feedback DAC.

**Noise Model:**

```matlab
% IDAC current noise
% File: sim_sense_noise.m, Line 32

[~, ref.idac.noise.idac_noise] = colouredNoise(...
    (obj.ctrl.N-1)*1/obj.imu.as.gyr.config.inf.fs, ...
    obj.imu.as.gyr.config.inf.fs, ...
    ref.idac.noise.freq, ...
    obj.ctrl.noise_sel('IDAC') * ref.idac.noise.ind, ...
    64345, 0, 0, 1, false);

% Total IDAC noise (including bandgap contribution)
idac_bg_noise = ref.idac.noise.idac_noise;
```

**IDAC to Rate Conversion:**

The IDAC noise appears at the rate integrator output through the feedback resistor:

$$V_{noise,rate} = I_{DAC,noise} \times R_{FB} \times \frac{K_{ref,noise}}{K_{in,noise}}$$

```matlab
% Bandgap+IDAC noise at C/V out for Rate integrator
C_sense.cvout_bg_noise = ...
    (obj.ctrl.noise_sel('Vref_SD') * idac_bg_noise * ...
     (sense.rate.coeff.Kref_noise / sense.rate.coeff.Kin_noise) * ...
     sense.rate.common.res_sig12) * ...
    obj.sys.rate_y.Value / obj.imu.as.gyr.config.inf.FSdps;
```

**Physical Insight:** The noise coefficient ratio $K_{ref,noise}/K_{in,noise}$ accounts for different integration of reference and signal through the pseudo-sinusoidal demodulation.

#### 3.2.3 Common-Mode Voltage LDO

The CM voltage LDO provides the virtual ground for the MEMS structure.

**Noise Components:**

```matlab
% CM voltage LDO noise
% File: sim_sense_noise.m, Line 42

[~, ref.cm.noise.ldo_noise] = colouredNoise(...
    (obj.ctrl.N-1)*1/obj.imu.as.gyr.config.inf.fs, ...
    obj.imu.as.gyr.config.inf.fs, ...
    ref.cm.noise.freq, ...
    obj.ctrl.noise_sel('VCM') * ref.cm.noise.vnd, ...
    87646, 0, 0, 1, false);
```

#### 3.2.4 Bandgap Coupling to CM Voltage

**Critical Correlation:** The CM LDO derives from the bandgap, creating correlated noise:

$$V_{CM,bg\_corr} = \frac{V_{CM}}{V_{BG}} \times V_{BG,noise}$$

```matlab
% Bandgap coupling to CM voltage
% File: sim_sense_noise.m, Line 51

ref.cm.noise.vref_noise = ...
    obj.ctrl.noise_sel('Vref_CM') * ...
    (ref.cm.vout/ref.bg.vout) * ref.bg.noise.bg_noise;
```

#### 3.2.5 Total CM Voltage Noise with Filtering

The CM voltage includes filtered LDO noise plus the bandgap-coupled component:

```matlab
% Combined CM noise with filtering
% File: sim_sense_noise.m, Line 54

[ref.cm.filt.b, ref.cm.filt.a] = sos2tf(ref.cm.filt.SOS, ref.cm.filt.G);

ref.cm.noise.cm_noise = ...
    filtfilt(ref.cm.filt.b, ref.cm.filt.a, ...
             ref.cm.noise.ldo_noise + ref.cm.noise.vref_noise) + ...
    (ref.cm.noise.thermal_vnd25k * sqrt(obj.imu.as.gyr.config.inf.fs/2) * ...
     randn(1,obj.ctrl.N));
```

**Physical Insight:** The `filtfilt()` function provides zero-phase filtering (typically ~2 kHz LPF), while thermal noise is added unfiltered.

#### 3.2.6 Uncorrelated CM Noise Component

For proper noise accounting, we must separate the correlated bandgap component:

```matlab
% CM voltage noise removing correlated bandgap part
C_sense.cvout_cm_noise_uncorrelated = ref.cm.noise.cm_noise - ref.cm.noise.vref_noise;
S_sense.cvout_cm_noise_uncorrelated = ref.cm.noise.cm_noise - ref.cm.noise.vref_noise;
Z_sense.cvout_cm_noise_uncorrelated = ref.cm.noise.cm_noise - ref.cm.noise.vref_noise;
```

**Why This Matters:** The bandgap noise appears through multiple paths (IDAC, CM voltage, potentially others). Separating the correlated component prevents double-counting when summing noise sources.

### 3.3 MEMS Noise Sources

#### 3.3.1 Series Resistor Noise

Similar to the drive channel, sense electrode traces contribute thermal noise.

**Calculation for Three Axes:**

```matlab
% MEMS sense series resistor noise
% File: sim_sense_noise.m, Line 59

sns_noise.cv.noise.C.vn_rs = ...
    obj.ctrl.noise_sel('R_MEMS_C') * ...
    sqrt(4*phys.kT*(obj.imu.gyro.gyro.par.Rs_GP1 + obj.imu.gyro.gyro.par.Rs_GN1));

sns_noise.cv.noise.S.vn_rs = ...
    obj.ctrl.noise_sel('R_MEMS_C') * ...
    sqrt(4*phys.kT*(obj.imu.gyro.gyro.par.Rs_GP2 + obj.imu.gyro.gyro.par.Rs_GN2));

sns_noise.cv.noise.Z.vn_rs = ...
    obj.ctrl.noise_sel('R_MEMS_C') * ...
    sqrt(4*phys.kT*(obj.imu.gyro.gyro.par.Rs_GP3 + obj.imu.gyro.gyro.par.Rs_GN3));
```

**Noise Gain with Distributed Parasitics:**

For sense electrodes, parasitic capacitances are distributed along the trace. The noise gain accounts for the average capacitance:

$$G_{n,rs} = \frac{(C_{G,CGM} + C_{G,SUB})/2}{C_{CVS}}$$

```matlab
% Noise gain (distributed parasitic model)
sns_noise.cv.noise.C.rs_gain = ...
    ((CG1_CGM + CG1_SUB)/2.0) / sense.cv.C.gyr_cap_cv_s;
```

#### 3.3.2 Brownian Mechanical Noise

Brownian motion of the sense mass creates a fundamental noise floor.

**Theory:**

The equivalent angular rate noise due to Brownian motion:

$$\Omega_{Brw} = \sqrt{\frac{4 k_B T \omega_0}{m_{eff} Q_{sense} S_{gain}^2 A_{drive}^2}}$$

where:
- $\omega_0$: Sense mode resonant frequency [rad/s]
- $m_{eff}$: Effective sense mass [kg]
- $Q_{sense}$: Sense mode quality factor
- $S_{gain}$: Coriolis gain (capacitance change per dps)
- $A_{drive}$: Drive amplitude [m]

**Implementation:**

```matlab
% MEMS Brownian noise
% File: sim_sense_noise.m, Line 74

C_sense.brw_noise_dps = ...
    obj.ctrl.noise_sel('MEMS_Brw') * ...
    obj.imu.gyro.gyro.fp.misc.brownian_noise.C * ...
    sqrt(obj.imu.as.gyr.config.inf.fs/2) * randn(1,obj.ctrl.N);
```

**Physical Insight:** Brownian noise is:
- White in frequency domain
- Proportional to √T
- Inversely proportional to √(m·Q)
- Reduced by higher drive amplitude
- A fundamental thermodynamic limit

**Typical Values:**
- High-performance gyroscope: 0.001-0.01 dps/√Hz
- Consumer gyroscope: 0.01-0.1 dps/√Hz
- Low-cost gyroscope: 0.1-1 dps/√Hz

### 3.4 Sense C/V Converter

#### 3.4.1 Sense C/V OTA Noise with Colored Spectrum

The sense C/V OTA is typically the dominant ASIC noise source. It includes thermal, flicker, and potentially RTN components.

**Total Capacitive Loading:**

$$C_{total} = C_{G,ALL} + C_{par,ASIC} + C_{in,pair} + C_{QC} + C_{CVS}$$

**Noise Gain:**

$$G_{n,CVS} = \frac{C_{total}}{C_{CVS}}$$

**Colored Noise Generation:**

```matlab
% Sense C/V OTA noise (colored)
% File: sim_sense_noise.m, Line 86

% Check if frequency-dependent noise data exists
if isfield(sense.cv.noise, 'vnd')
    cv_noise_vnd = sense.cv.noise.vnd;  % Frequency-dependent
else
    cv_noise_vnd = sense.cv.noise.Vn_opa * ones(1, length(sense.cv.noise.freq));
end

[~, sns_noise.cv.noise.cvs_noise] = colouredNoise(...
    (obj.ctrl.N-1)*1/(obj.imu.as.gyr.config.inf.fs), ...
    obj.imu.as.gyr.config.inf.fs, ...
    sense.cv.noise.freq, ...
    obj.ctrl.noise_sel('CVS') * cv_noise_vnd, ...
    542353, 0, 0, 1, false);
```

**Pseudo-Sinusoidal Modulation:**

The noise is explicitly modulated by the pseudo-sinusoidal demodulation waveform:

```matlab
% Generate pseudo-sin signal
pseudo_sin_signal = repmat(obj.pseudo_sin(0:ceil(obj.ctrl.N/192)*192-1, 192), ...
                           ceil(obj.ctrl.N/192));

% Apply pseudo-sin modulation to noise
C_sense.cvout_cvs_ota = ...
    (sns_noise.cv.noise.cvs_noise * sns_noise.cv.noise.C.cvs_gain) .* ...
    pseudo_sin_signal(1:obj.ctrl.N);
```

**Physical Insight:** The pseudo-sin modulation:
- Translates noise from baseband to drive frequency
- Creates sidebands around f_drive
- Only demodulated components contribute to output
- Effectively band-limits the noise

#### 3.4.2 RTN in Input Transistors

Random Telegraph Noise from oxide traps in input transistors can be significant in small-geometry MOSFETs.

**See Section 7 for complete RTN modeling details.**

**Quick Integration Example:**

```matlab
% Add RTN to sense C/V OTA (optional enhancement)
if obj.ctrl.noise_sel('CVS_RTN')
    % Generate RTN signal
    num_traps = 3;
    tau_high = [10e-6, 25e-6, 50e-6];  % Trap time constants
    tau_low = [8e-6, 20e-6, 45e-6];
    amplitude = [1e-3, 2e-3, 1.5e-3];  % mV steps
    
    rtn_signal = generate_rtn(obj.ctrl.N, obj.imu.as.gyr.config.inf.fs, ...
                              amplitude, tau_high, tau_low, num_traps, seed);
    
    % Add to OTA noise with appropriate gain
    C_sense.cvout_cvs_ota = C_sense.cvout_cvs_ota + ...
                            rtn_signal * sns_noise.cv.noise.C.cvs_gain;
end
```

### 3.5 Rate Integrator

#### 3.5.1 Signal Resistor Noise

Similar to the phase integrator, the rate integrator uses pseudo-sinusoidal switching.

**Time-Averaged Noise:**

```matlab
% Rate integrator signal resistor noise
% File: sim_sense_noise.m, Line 103

sns_noise.rate.noise.common.vn_rsig = ...
    obj.ctrl.noise_sel('Rint_Rsign') * ...
    sqrt((2*4*phys.kT*sense.rate.common.res_sig12*sense.rate.phi_hi / ...
          (sense.rate.phi_per/4)) + ...
         (2*4*phys.kT*sense.rate.common.res_sig1*sense.rate.phi_lo / ...
          (sense.rate.phi_per/4)));
```

#### 3.5.2 Rate Integrator OTA Noise

**Noise Scaling:**

```matlab
% Rate integrator OTA noise
% File: sim_sense_noise.m, Line 116

C_sense.cvout_rate_ota = ...
    (sqrt((obj.ctrl.noise_sel('Rint_opa')*sense.rate.noise.Vn_opa)^2 + ...
          sns_noise.rate.noise.C.vint_ota_gain^2) * ...
     sense.rate.coeff.Kin / sense.rate.coeff.Kin_noise) * ...
    sqrt(obj.imu.as.gyr.config.inf.fs/2) * randn(1,obj.ctrl.N);
```

#### 3.5.3 Finite OTA Gain Effects

Virtual ground error due to finite OTA gain can modulate the signal, creating signal-dependent noise.

**Error Voltage:**

$$V_{\epsilon} = \frac{V_{quad}}{10^{A_{OL}/20}}$$

where $A_{OL}$ is the open-loop gain at the signal frequency (typically 70 kHz).

**Resulting Noise:**

```matlab
% Finite OTA gain error
sns_noise.rate.noise.C.veps = ...
    qc.C.quad_res * (sense.cv.C.cv_out/obj.imu.as.gyr.config.inf.FSdps) / ...
    10^(sense.rate.common.gain_70kHz/20);

% Noise contribution (currently disabled in model)
sns_noise.rate.noise.C.vint_ota_gain = ...
    0 * obj.ctrl.noise_sel('Rint_opa_gain') * ...
    sns_noise.rate.noise.C.veps * sense.rate.C.FS * ...
    (2*sqrt((1-(2/pi))*2/obj.imu.as.gyr.config.inf.fs)) / (pi*obj.imu.as.gyr.cm.vref);
```

**Note:** This term is currently disabled (factor of 0) because for typical OTA gains >80 dB, it's negligible compared to other sources. Enable if OTA gain is marginal.

### 3.6 ASIC Parallel Resistor Noise

Similar to drive channel, large parallel resistors contribute noise.

**Impedance at Sense Frequency:**

The impedance is calculated at the *drive frequency* where the signal is modulated:

$$Z_{in} = \frac{1}{j\omega_{drive} C_{total}}$$

**Current to Voltage Conversion:**

```matlab
% ASIC parallel resistor noise
% File: sim_sense_noise.m, Line 136

sns_noise.cv.noise.C.vn_rp = ...
    obj.ctrl.noise_sel('Rpar_C') * ...
    sqrt(2*4*phys.kT/obj.imu.as.gyr.par.Rpar_in) * ...
    sense.cv.C.cv_gain_Iin_Vout;

C_sense.cvout_rp = ...
    sns_noise.cv.noise.C.vn_rp * ...
    sqrt(obj.imu.as.gyr.config.inf.fs/2) * randn(1,obj.ctrl.N);
```

---

## 4. Common-Mode Voltage Noise Propagation

### 4.1 Overview

Despite the fully differential architecture, common-mode voltage noise couples into the signal path through several mechanisms. These coupling paths can be significant contributors to total noise.

**Four Primary Coupling Mechanisms:**

1. **CM noise over C₀ offset**: Capacitor mismatch breaks differential symmetry
2. **CM noise through parasitic mismatch**: Substrate parasitics differ between electrodes
3. **CM noise via rate signal modulation**: Signal-dependent coupling
4. **CM noise demodulated by quadrature**: Residual quadrature provides carrier

### 4.2 CM Noise Over C₀ Capacitor Offset

#### 4.2.1 Physical Mechanism

The differential sense capacitors are never perfectly matched. An offset $\Delta C_0 = C_{GP} - C_{GN}$ exists due to:
- Fabrication tolerances (±0.1-1% typical)
- Comb finger geometry variations
- Drive amplitude asymmetry
- Mechanical stress

When common-mode voltage varies, this offset converts CM voltage into differential voltage:

$$V_{diff} = V_{CM} \times \frac{\Delta C_0}{C_{total}}$$

#### 4.2.2 Noise Transfer to C/V Output

**Sinusoidal Modulation:**

The drive motion amplitude-modulates the offset, creating a signal at the drive frequency:

$$e_{n,CM,C0} = V_{CM,noise} \times \frac{\Delta C_0}{C_{CVS}} \times \frac{1}{\sqrt{2}} \times \sin(2\pi f_d t)$$

**Key Points:**
- Modulated at $f_d$ because offset scales with drive amplitude
- Factor $1/\sqrt{2}$: Only in-phase component demodulates into rate
- Demodulation occurs in rate integrator

**Implementation:**

```matlab
% CM voltage noise over C0 capacitor offset
% File: sim_sense_noise.m, Line 128

C_sense.cvout_cm_c0_noise = ...
    obj.ctrl.noise_sel('CM_C0') * ...
    (ref.cm.noise.cm_noise * ...
     (obj.imu.gyro.gyro.par.CGP1_CGN1/sense.cv.C.gyr_cap_cv_s)) / sqrt(2) .* ...
    sin(2*pi*obj.imu.gyro.gyro.stat.mc.fd*obj.ctrl.t);

% Similar for S and Z axes
S_sense.cvout_cm_c0_noise = ...
    obj.ctrl.noise_sel('CM_C0') * ...
    (ref.cm.noise.cm_noise * ...
     (obj.imu.gyro.gyro.par.CGP2_CGN2/sense.cv.S.gyr_cap_cv_s)) / sqrt(2) .* ...
    sin(2*pi*obj.imu.gyro.gyro.stat.mc.fd*obj.ctrl.t);
```

#### 4.2.3 Design Implications

**Sensitivity:**

For typical values:
- $\Delta C_0 / C_{nom}$ = 0.5% (5000 ppm)
- $V_{CM,noise}$ = 100 nV/√Hz
- $C_{CVS}$ = 1 pF

Equivalent input noise ≈ 5 nV/√Hz

**Mitigation Strategies:**
1. Tight capacitor matching (layout symmetry)
2. Differential drive excitation (balances offsets)
3. Trim/calibration of C₀ offset
4. Higher feedback capacitor $C_{CVS}$ (reduces gain)

### 4.3 CM Noise Through Parasitic Capacitance Mismatch

#### 4.3.1 Physical Mechanism

The MEMS proof mass connects to the substrate through parasitic capacitances. Mismatch in these parasitics ($C_{GP,SUB} \neq C_{GN,SUB}$) creates an additional CM coupling path.

**Equivalent Circuit:**

```
        VCM
         |
    CGP,SUB  CGN,SUB
         |    |
         GP   GN  (Sense electrodes)
         |    |
        CVS  CVS
         |____|
           |
          Vout
```

#### 4.3.2 Noise Transfer Function

**Parasitic Divider Effect:**

The mismatch creates a CM-to-differential conversion:

$$G_{CM,par} = \frac{C_{G,CGM}}{C_{G,CGM} + C_{G,SUB} + C_{CVS} + \Delta C_{SUB}/2}$$

**Total Noise at C/V Output:**

$$e_{n,CM,par} = G_{CM,par} \times \frac{\Delta C_{SUB}}{C_{CVS}} \times V_{CM,uncorr} \times \frac{1}{\sqrt{2}} \times \sin(2\pi f_d t)$$

where:
- $\Delta C_{SUB} = C_{GP,SUB} - C_{GN,SUB}$: Parasitic mismatch
- $V_{CM,uncorr}$: CM noise excluding correlated bandgap component

**Implementation:**

```matlab
% CM parasitic mismatch coupling gain
% File: sim_sense_noise.m, Line 139

sns_noise.cv.noise.C.cm_cpar_gain = ...
    CG1_CGM / (CG1_CGM + CG1_SUB + sense.cv.C.gyr_cap_cv_s + dCG1_SUB/2);

% CM noise through parasitic mismatch
C_sense.cvout_cm_cpar_noise = ...
    obj.ctrl.noise_sel('CM_Par') * ...
    sns_noise.cv.noise.C.cm_cpar_gain * (dCG1_SUB)/sense.cv.C.gyr_cap_cv_s * ...
    C_sense.cvout_cm_noise_uncorrelated / sqrt(2) .* ...
    sin(2*pi*obj.imu.gyro.gyro.stat.mc.fd*obj.ctrl.t);
```

#### 4.3.3 Why Use Uncorrelated CM Noise?

We use `cvout_cm_noise_uncorrelated` because:
1. The bandgap-correlated component already appears through IDAC path
2. Using total CM noise would double-count the bandgap contribution
3. Only the LDO-generated noise couples through this path

#### 4.3.4 Design Implications

**Sensitivity:**

For typical values:
- $\Delta C_{SUB}$ = 1-10 fF (mismatch)
- $C_{G,SUB}$ = 100-500 fF (total)
- $V_{CM,noise}$ = 100 nV/√Hz

Requires **CMRR > 80 dB** at drive frequency.

**Mitigation Strategies:**
1. Symmetric substrate contact layout
2. Guard rings to control parasitics
3. Differential substrate biasing
4. Shield electrodes to define parasitic paths

### 4.4 CM Noise Through Rate Signal Modulation

#### 4.4.1 Physical Mechanism

The Coriolis-induced capacitance change is:

$$\Delta C = S_{dd} \times \Omega_{rate} \times \sin(2\pi f_d t)$$

This capacitance variation causes the *voltage* across it to modulate the effective capacitance seen by the C/V converter.

#### 4.4.2 Signal-Dependent Noise Coupling

**Voltage-Dependent Capacitance:**

When CM voltage varies, the capacitance change translates to a charge variation:

$$\Delta Q = \Delta C \times V_{CM,noise} = S_{dd} \times \Omega_{rate} \times V_{CM,noise}$$

**At C/V Output:**

$$e_{n,CM,rate} = V_{CM,uncorr} \times \frac{S_{dd} \times \Omega_{rate}}{C_{CVS}}$$

**Implementation:**

```matlab
% CM voltage noise through rate signal
% File: sim_sense_noise.m, Line 166

C_sense.cvout_cm_rate_noise = ...
    obj.ctrl.noise_sel('CM_rate') * ...
    (C_sense.cvout_cm_noise_uncorrelated * ...
     (obj.imu.gyro.gyro.stat.mc.Sdd_C_Fdps * obj.sys.rate_x.Value / ...
      sense.cv.C.gyr_cap_cv_s));
```

#### 4.4.3 Key Observations

**Signal-Dependent:**
- Proportional to applied rate
- Zero at zero rate
- Can be positive or negative depending on rate direction

**Not Modulated:**
- Already at baseband (DC signal)
- Directly appears in rate integrator
- No $\sin(2\pi f_d t)$ term needed

**Implications for Testing:**
- Difficult to separate from other noise in static tests
- Changes with applied rate
- Requires CM voltage noise characterization

### 4.5 CM Noise Demodulated by Quadrature Motion

#### 4.5.1 Physical Mechanism

Residual quadrature motion (after electronic compensation) creates an alternative carrier for demodulating CM noise.

**Quadrature Motion:**

$$x_{quad}(t) = Q_{res} \times A_{drive} \times \sin(2\pi f_d t + \phi_{quad})$$

where $Q_{res}$ is the residual quadrature in dps-equivalent after compensation.

**Capacitance Change:**

$$\Delta C_{quad} = S_{dd} \times Q_{res} \times \sin(2\pi f_d t + \phi_{quad})$$

#### 4.5.2 Demodulation of CM Noise

The mechanical phase angle $\phi_{mech}$ determines how much CM noise appears in rate vs. quadrature channels:

**Rate Component:**

$$e_{n,CM,quad,rate} = \frac{S_{dd} \times Q_{res}}{C_{CVS}} \times \sin(\phi_{mech}) \times V_{CM,noise}$$

**Quadrature Component:**

$$e_{n,CM,quad,quad} = \frac{S_{dd} \times Q_{res}}{C_{CVS}} \times \cos(\phi_{mech}) \times V_{CM,noise}$$

**Implementation:**

```matlab
% CM noise demodulated by quadrature - Rate component
% File: sim_sense_noise.m, Line 150

C_sense.cvout_cm_qdem_noise_r = ...
    obj.ctrl.noise_sel('CM_Qdem') * ...
    ((obj.imu.gyro.gyro.stat.mc.Sdd_C_Fdps * ...
      (obj.imu.gyro.gyro.stat.mc.quad_C + qc.C.quad_res)) / ...
     sense.cv.C.gyr_cap_cv_s) * ...
    sin(obj.imu.gyro.gyro.fp.detection.mechanical_phase.C/1000) * ...
    ref.cm.noise.cm_noise;

% Quadrature component
C_sense.cvout_cm_qdem_noise_q = ...
    obj.ctrl.noise_sel('CM_Qdem') * ...
    ((obj.imu.gyro.gyro.stat.mc.Sdd_C_Fdps * ...
      (obj.imu.gyro.gyro.stat.mc.quad_C + qc.C.quad_res)) / ...
     sense.cv.C.gyr_cap_cv_s) * ...
    cos(obj.imu.gyro.gyro.fp.detection.mechanical_phase.C/1000) * ...
    ref.cm.noise.cm_noise;
```

#### 4.5.3 Design Implications

**Key Points:**
- Proportional to residual quadrature (after electronic compensation)
- Mechanical phase determines rate/quad split
- Uses *total* CM noise (not uncorrelated) because it's a direct modulation

**Typical Values:**
- $Q_{res}$ after compensation: 1-50 dps equivalent
- $\phi_{mech}$: 85-95° (close to quadrature)
- Therefore: mostly appears in rate channel

**Mitigation:**
1. Better quadrature compensation (reduce $Q_{res}$)
2. Mechanical trimming of quadrature
3. Lower CM voltage noise
4. Higher $C_{CVS}$ (reduces gain)

### 4.6 Quadrature Compensation Coupling

#### 4.6.1 Mechanism

The quadrature compensation (QC) injects a signal from the drive C/V into the sense C/V through a coupling capacitor $C_{QC}$.

**Noise Feedthrough:**

Drive C/V noise couples directly:

$$e_{n,QC} = e_{n,drive,CVD} \times \frac{C_{QC}}{C_{CVS}}$$

**Implementation:**

```matlab
% Noise from drive C/V through quadrature compensation
% File: sim_sense_noise.m, Line 174

C_sense.cvout_noise_qc = ...
    obj.ctrl.noise_sel('Quad_Comp') * ...
    drive.noise.cv.noise.vn_rand * qc.C.gyr_cap_qc / sense.cv.C.gyr_cap_cv_s;
```

#### 4.6.2 Design Trade-offs

**QC Capacitor Sizing:**
- Larger $C_{QC}$: More compensation range, but more noise coupling
- Smaller $C_{QC}$: Less noise, but limited compensation range

**Typical Values:**
- $C_{QC}/C_{CVS}$ = 0.01-0.1
- Noise coupling = -20 to -40 dB

### 4.7 Summary of CM Coupling Paths

| Coupling Path | Modulated? | Signal-Dependent? | Typical Contribution |
|---------------|------------|-------------------|---------------------|
| C₀ offset | Yes (sin ωdt) | No | 5-20% |
| Parasitic mismatch | Yes (sin ωdt) | No | 2-10% |
| Rate modulation | No | Yes (∝ rate) | <5% at ±100 dps |
| Quad demodulation | No | Yes (∝ Qres) | 1-5% |
| QC coupling | Yes (sin ωdt) | No | 2-15% |

**Total CM Contribution:** Can be 15-40% of total noise if not properly managed.

---
