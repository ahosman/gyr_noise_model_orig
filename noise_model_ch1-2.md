# Open-Loop Capacitive MEMS Gyroscope Noise Model
## Complete Documentation with Implementation Details

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
