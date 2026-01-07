# Open-Loop Capacitive MEMS Gyroscope Noise Model
## Complete Documentation with Implementation Details

**Document Version:** 3.0  
**Last Updated:** December 2024  
**Author:** Ahmed Hussein (BST/EIC6)  
**MATLAB Implementation:** BAI482 Noise Model

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Drive Channel Noise](#2-drive-channel-noise)
3. [Sense Channel Noise](#3-sense-channel-noise)
4. [Common-Mode Voltage Noise Propagation](#4-common-mode-voltage-noise-propagation)
5. [Sigma-Delta Modulator](#5-sigma-delta-modulator)
6. [Pseudo-Sinusoidal Demodulation](#6-pseudo-sinusoidal-demodulation)
7. [Random Telegraph Noise (RTN)](#7-random-telegraph-noise-rtn)
8. [Noise Budgeting and Analysis](#8-noise-budgeting-and-analysis)
9. [Implementation Guide](#9-implementation-guide)
10. [Validation and Testing](#10-validation-and-testing)
11. [References](#11-references)
12. [Appendices](#12-appendices)

---

## 1. Introduction

### 1.1 System Overview

This document describes the comprehensive noise model for an open-loop capacitive MEMS gyroscope system. The model encompasses all significant noise sources from the MEMS structure through the ASIC readout chain to the digital output.

**Key Features:**
- Three-axis gyroscope (C, S, Z axes corresponding to Roll, Pitch, Yaw)
- Capacitive sense with pseudo-sinusoidal demodulation
- MASH 1-1 Sigma-Delta modulator architecture
- Comprehensive common-mode noise coupling analysis
- Random Telegraph Noise (RTN) modeling
- Temperature-dependent noise sources

### 1.2 Noise Model Scope

The model includes:

**Physical Noise Sources:**
- Thermal (Johnson-Nyquist) noise from resistors
- Brownian mechanical noise from MEMS structure
- Flicker (1/f) noise from active devices
- Random Telegraph Noise (RTN) from oxide traps
- Shot noise from current sources
- Phase noise from VCO and PLL

**Coupling Mechanisms:**
- Direct signal path noise
- Common-mode voltage coupling (4 distinct paths)
- Quadrature compensation feedthrough
- Phase noise modulation of sigma-delta feedback
- Parasitic capacitance mismatch effects

**Signal Processing:**
- Capacitance-to-voltage conversion
- Pseudo-sinusoidal demodulation
- Rate integration with feedback
- Two-stage sigma-delta modulation
- Digital filtering and decimation

### 1.3 Coordinate Systems and Conventions

```
MEMS Gyroscope Coordinate System:
  C-axis (Channel 1): Roll  - Rotation about lateral axis
  S-axis (Channel 2): Pitch - Rotation about longitudinal axis  
  Z-axis (Channel 3): Yaw   - Rotation about vertical axis

Sign Convention:
  - Positive rotation: Right-hand rule
  - Capacitance increase: Positive signal
  - Drive motion: Sinusoidal at f_d ≈ 35 kHz (typ)
```

### 1.4 File Organization

```matlab
% Main class definition
classdef noise_model < handle
    properties
        ctrl % Simulation controls
        imu  % IMU/gyroscope object
        sys  % System state
        info % Results structure
        spectrum % Spectral analysis
        adev % Allan deviation
    end
    
    methods
        function obj = noise_model(imu, sys, noise_sel)
            % Constructor: Initialize and run noise simulation
            obj.imu  = imu;
            obj.sys  = sys;
            obj.ctrl = noise_model.init_ctrl();
            obj.ctrl.noise_sel = noise_model.init_noise_cfg(noise_sel);
            
            % Run noise simulations
            obj.info.drive.noise = obj.sim_drive_noise;
            [obj.info.sense.noise, ...
             obj.info.C.sense.noise, ...
             obj.info.S.sense.noise, ...
             obj.info.Z.sense.noise, ...
             obj.spectrum.C.sense, ...
             obj.spectrum.S.sense, ...
             obj.spectrum.Z.sense] = obj.sim_sense_noise;
            [obj.spectrum.C.rate, ...
             obj.spectrum.S.rate, ...
             obj.spectrum.Z.rate] = obj.sim_rate_sd_noise;
        end
    end
end
```

---

## 2. Drive Channel Noise

### 2.1 Drive Channel Architecture

The drive channel maintains the MEMS resonator oscillation at its resonant frequency (typically 35 kHz). It consists of:

1. **Drive C/V Converter:** Converts drive electrode capacitance to voltage
2. **Phase-Locked Loop (PLL):** Maintains oscillation phase and amplitude
3. **Phase Integrator:** Demodulates drive motion for phase control
4. **VCO:** Generates drive excitation signal

```
MEMS Drive    Drive C/V      Phase         VCO         Drive
Electrodes → Converter  →  Integrator  → Controller → Excitation
   (DP/DN)     (CVD)         (Pint)                      (Vdrive)
      ↑______________________________________________|
                    Closed-Loop Drive Control
```

### 2.2 Drive C/V Converter Noise

#### 2.2.1 OTA Noise

The drive C/V operational transconductance amplifier (OTA) contributes thermal and flicker noise.

**Noise Density at C/V Output:**

$$e_{n,out,cvd} = e_{n,opa} \times \frac{C_{0,DRV} + C_{p,DRV,MEMS} + C_{p,DRV,ASIC} + C_{CVD}}{C_{CVD}}$$

where:
- $e_{n,opa}$: Input-referred noise of Drive C/V amplifier [V/√Hz]
- $C_{0,DRV}$: Active drive capacitance (DP/DN to CM) [F]
- $C_{p,DRV,MEMS}$: MEMS parasitic capacitance on DP/DN [F]
- $C_{p,DRV,ASIC}$: ASIC parasitic capacitance on DP/DN [F]
- $C_{CVD}$: Drive C/V integration capacitance [F]

**Implementation:**

```matlab
% Drive C/V OTA noise calculation
% File: sim_drive_noise.m, Line 19

C0_D = (obj.imu.gyro.gyro.stat.mc.Cdp + obj.imu.gyro.gyro.stat.mc.Cdn)/2.0;

drv_noise.cv.noise.vn_cvd = ...
    obj.ctrl.noise_sel('CVD') * drive.cv.noise.Vn_opa * ...
    (C0_D + (obj.imu.gyro.gyro.par.Cpar_DP + obj.imu.gyro.gyro.par.Cpar_DN)/2 + ...
     obj.imu.as.gyr.par.Cpar_C + drive.cv.gyr_cap_cv_d) / drive.cv.gyr_cap_cv_d;
```

**Typical Values:**
- $e_{n,opa}$ ≈ 5-10 nV/√Hz
- Noise gain ≈ 5-20 (depending on parasitic loading)
- Output noise ≈ 25-200 nV/√Hz

#### 2.2.2 MEMS Series Resistor Noise

The aluminum or gold traces connecting the drive electrodes contribute thermal noise.

**Noise Density (Two Resistors in Series):**

$$e_{n,out,rs} = \sqrt{2 \times 4 k_B T \times Rs_D} \times \frac{C_{0,DRV} + C_{p,DRV,MEMS}}{C_{CVD}}$$

where:
- $k_B$ = 1.38×10⁻²³ J/K: Boltzmann's constant
- $T$: Absolute temperature [K]
- $Rs_D$: Resistance of one drive electrode [Ω]

**Physical Insight:** Only the capacitance *before* the series resistor contributes to the noise gain. The feedback capacitor $C_{CVD}$ appears in the denominator as it forms the integrator.

**Implementation:**

```matlab
% MEMS series resistor noise
% File: sim_drive_noise.m, Line 26

% Thermal noise density of series resistors
drv_noise.cv.noise.vn_rsd = ...
    obj.ctrl.noise_sel('R_MEMS_D') * ...
    sqrt(4*phys.kT*(obj.imu.gyro.gyro.par.Rs_DP + obj.imu.gyro.gyro.par.Rs_DN));

% Noise gain (only C before Rs counts)
drv_noise.cv.noise.rs_gain = ...
    (C0_D + (obj.imu.gyro.gyro.par.Cpar_DP + obj.imu.gyro.gyro.par.Cpar_DN)/2) / ...
    drive.cv.gyr_cap_cv_d;

% Total noise at C/V output
drv_noise.cv.noise.vn_rs = ...
    drv_noise.cv.noise.rs_gain * drv_noise.cv.noise.vn_rsd;
```

**Typical Values:**
- $Rs_D$ ≈ 50-200 Ω per electrode
- Thermal noise ≈ 1-2 nV/√Hz
- With gain ≈ 5-50 nV/√Hz at output

#### 2.2.3 ASIC Parallel Resistor Noise

Large resistors (1-10 MΩ) provide DC bias and ESD protection. Their noise appears as current noise.

**Noise Density at C/V Output:**

$$e_{n,out,rp} = \sqrt{\frac{2 \times 4 k_B T}{Rp_{in}}} \times \frac{1}{2\pi f_{drv} \times C_{CVD}}$$

where:
- $Rp_{in}$: Parallel input resistor (one of two) [Ω]
- $f_{drv}$: Drive frequency [Hz]
- Factor of 2: Two resistors (DP and DN)

**Physical Insight:** The parallel resistor noise appears as current noise. The impedance at the drive frequency determines the voltage noise: $Z = 1/(j\omega C_{CVD})$.

**Implementation:**

```matlab
% ASIC parallel resistor noise
% File: sim_drive_noise.m, Line 40

% Current noise spectral density: sqrt(4kT/R) [A/√Hz]
% Impedance at drive frequency: 1/(2πf·C) [Ω]
drv_noise.cv.noise.vn_rp = ...
    obj.ctrl.noise_sel('Rpar_D') * ...
    sqrt(2*4*phys.kT/obj.imu.as.gyr.par.Rpar_in) / ...
    (2*pi*obj.imu.gyro.gyro.stat.mc.fd*drive.cv.gyr_cap_cv_d);
```

**Typical Values:**
- $Rp_{in}$ ≈ 1-10 MΩ
- At 35 kHz with 1 pF: ≈ 10-30 nV/√Hz

#### 2.2.4 Total Drive C/V Noise

**RMS Combination:**

$$e_{n,drive,cvd} = \sqrt{e_{n,out,cvd}^2 + e_{n,out,rs}^2 + e_{n,out,rp}^2}$$

**Implementation:**

```matlab
% Total noise density of Drive C/V at output
drv_noise.cv.noise.vn = ...
    sqrt(drv_noise.cv.noise.vn_cvd^2 + ...
         drv_noise.cv.noise.vn_rs^2 + ...
         drv_noise.cv.noise.vn_rp^2);

% Generate noise time series
drv_noise.cv.noise.vn_rand = ...
    drv_noise.cv.noise.vn * sqrt(obj.imu.as.gyr.config.inf.fs/2) * randn(1,obj.ctrl.N);
```

### 2.3 Phase-Locked Loop Noise

#### 2.3.1 Phase Integrator Resistor Noise

The phase integrator uses pseudo-sinusoidal demodulation with switched resistors.

**Effective Resistance:**

The pseudo-sinusoidal demodulation uses two resistor values ($R_{sig,1} = 7R_{unit}$ and $R_{sig,2} = 5R_{unit}$) switched at different phases:

$$R_{sig,12} = R_{sig,1} \parallel R_{sig,2} = \frac{7R_{unit} \times 5R_{unit}}{7R_{unit} + 5R_{unit}} = \frac{35}{12}R_{unit}$$

**Time-Averaged Noise:**

$$e_{n,pint,rsig} = \sqrt{2 \times 4 k_B T \times \left[\frac{\phi_{hi}R_{sig,12} + \phi_{lo}R_{sig,1}}{\phi_{per}/4}\right]}$$

where:
- $\phi_{hi}$: High phase duration (both resistors on)
- $\phi_{lo}$: Low phase duration (only $R_{sig,1}$ on)
- $\phi_{per}$: Total pseudo-sin period (192 samples typical)
- Division by $\phi_{per}/4$: Normalizes to quarter period

**Implementation:**

```matlab
% Phase integrator signal resistor noise
% File: sim_drive_noise.m, Line 51

drv_noise.pll.noise.vn_rsig = ...
    obj.ctrl.noise_sel('Pint_Rsign') * ...
    sqrt((2*4*phys.kT*drive.pll.res_sig12*drive.pll.phi_hi/(drive.pll.phi_per/4)) + ...
         (2*4*phys.kT*drive.pll.res_sig*drive.pll.phi_lo/(drive.pll.phi_per/4)));
```

**Typical Values:**
- $R_{unit}$ ≈ 50-200 kΩ
- $\phi_{hi}/\phi_{per}$ ≈ 0.25-0.5
- Effective noise ≈ 10-40 nV/√Hz

#### 2.3.2 Phase Integrator OTA Noise

The phase integrator OTA noise is scaled by the charge integration coefficient.

**Charge Integration:**

The integrated charge depends on the pseudo-sinusoidal correlation:

$$Q_{1A} = \frac{A_1 I_{r0}}{2\pi f_d}\left[-\cos\left(2\pi\frac{\phi_{lo}+\phi_{grd}}{\phi_{per}}\right) + \cos\left(2\pi\frac{\phi_{grd}}{\phi_{per}}\right)\right]$$

$$Q_{1B} = \frac{I_{r0}}{2\pi f_d}\left[-\cos\left(2\pi\frac{\phi_{hi}+\phi_{lo}+\phi_{grd}}{\phi_{per}}\right) + \cos\left(2\pi\frac{\phi_{lo}+\phi_{grd}}{\phi_{per}}\right)\right]$$

**Signal Gain Coefficient:**

$$K_{in} = 4 \times [Q_{1A} + Q_{1B}]$$

**Noise Scaling:**

$$e_{n,pint,opa} = e_{n,opa} \times \frac{K_{in}}{K_{n,in}}$$

where $K_{n,in} \approx 0.7613$ is the noise coefficient extracted from FFT analysis.

**Implementation:**

```matlab
% Phase integrator OTA noise
% File: sim_drive_noise.m, Line 64

drv_noise.pll.noise.vn_int = ...
    obj.ctrl.noise_sel('Pint_opa') * drive.pll.noise.Vn_opa * ...
    drive.pll.coeff.Kin / drive.pll.coeff.Kin_noise;
```

**Physical Insight:** The signal coefficient $K_{in}$ differs from the noise coefficient $K_{n,in}$ because:
- Signal: Coherent integration over full period
- Noise: RMS integration with different weighting

This is a fundamental property of pseudo-sinusoidal demodulation (see Section 6 for details).

#### 2.3.3 VCO Phase Noise

The VCO generates colored phase noise with a characteristic $1/f^3$ profile close to the carrier.

**Phase Noise Specification:**

Typically specified at discrete frequency offsets from the carrier:

| Offset from Carrier | Phase Noise [dBc/Hz] |
|---------------------|----------------------|
| 100 Hz              | -60 to -80           |
| 1 kHz               | -80 to -100          |
| 10 kHz              | -100 to -120         |
| 100 kHz             | -120 to -140         |

**Colored Noise Generation:**

```matlab
% VCO phase noise (colored)
% File: sim_drive_noise.m, Line 75

[~, phi_noise] = colouredNoise(...
    (obj.ctrl.N-1)/obj.imu.as.gyr.config.inf.fs, ...  % Duration [s]
    obj.imu.as.gyr.config.inf.fs, ...                  % Sampling freq [Hz]
    drive.pll.noise.pll_freq, ...                      % Frequency points [Hz]
    (10.^(drive.pll.noise.pll_phnoise/20)) * ...     % PSD values [rad/√Hz]
        drive.pll.noise.scale_noise/4/sqrt(2) * sqrt(phys.T/300), ...
    12574, 0, 0, 0, false);                           % Seed and options

drv_noise.pll.noise.phase_noise_coloured = ...
    obj.ctrl.noise_sel('phase_noise_VCO') * phi_noise;
```

**Temperature Scaling:** Phase noise typically improves with temperature as $\sqrt{T_{ref}/T}$ due to reduced thermal agitation.

#### 2.3.4 Phase Noise Front-End Filtering

The drive C/V noise contributes to phase noise through the phase detector. This is filtered by a low-pass filter before combining with VCO noise.

**Implementation:**

```matlab
% Front-end phase noise filtering
% File: sim_drive_noise.m, Line 89

Htest = load('filt_200hz.mat');  % 200 Hz LPF

% Filter the white noise contribution
drv_noise.pll.noise.phase_noise_fe_v = ...
    filter(Htest.Hlp, ...
           drv_noise.pll.noise.phase_noise_coloured * ...
           sqrt(obj.imu.as.gyr.config.inf.fs/2) .* ...
           drv_noise.cv.noise.vn_rand);
```

**Physical Insight:** The PLL loop filter bandwidth (typically 100-500 Hz) limits how fast the phase can track. High-frequency noise from the drive C/V is attenuated by this filtering.

#### 2.3.5 Total Drive Phase Noise

**Combined Phase Noise:**

```matlab
% Total drive phase noise
% Optional: Add discrete tone at drive frequency
drv_noise.pll.peak_at_25k = 0*10^(-70/20)/64;  % Typically disabled

drv_noise.pll.noise.phase_noise_drive = ...
    sqrt(2) * ((drv_noise.pll.peak_at_25k * ...
                sin(2*pi*obj.imu.gyro.gyro.stat.drive.tuned_frequency.mean*obj.ctrl.t) + ...
                drv_noise.pll.noise.phase_noise_fe_v + ...
                drv_noise.pll.noise.phase_noise_coloured));
```

**Conversion to Phase Noise [dBc]:**

$$\phi_{noise}[dBc/Hz] = 20 \log_{10}\left(\frac{e_{n,pll}}{K_{in} \times V_{cv,out} \times 4\sqrt{2}}\right)$$

```matlab
drv_noise.pll.noise.phase_noise = ...
    drv_noise.pll.noise.vn_pll / (drive.pll.coeff.Kin*drive.cv.vout) / 4 / sqrt(2);

drv_noise.pll.noise.phase_noise_dBc = 20*log10(drv_noise.pll.noise.phase_noise);
```

### 2.4 Drive Noise Visualization

**Optional Spectral Plots:**

```matlab
if obj.ctrl.plot_spectra_d == 1
    figure(101);
    pwelch(drv_noise.pll.noise.phase_noise_drive,[],[],[],obj.imu.as.gyr.config.inf.fs);
    pwelch(drv_noise.pll.noise.phase_noise_fe_v,[],[],[],obj.imu.as.gyr.config.inf.fs);
    pwelch(drv_noise.pll.noise.phase_noise_coloured,[],[],[],obj.imu.as.gyr.config.inf.fs);
    set(gca,'xscale','log');
    legend('Phase Noise Drive', 'Phase Noise FE V', 'Phase Noise Coloured');
end
```

---

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
