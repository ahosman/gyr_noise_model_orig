
## 11. References

### 11.1 MEMS Gyroscope Theory

1. Yazdi, N., Ayazi, F., & Najafi, K. (1998). "Micromachined inertial sensors." *Proceedings of the IEEE*, 86(8), 1640-1659.

2. Acar, C., & Shkel, A. (2008). *MEMS Vibratory Gyroscopes: Structural Approaches to Improve Robustness*. Springer.

3. Geen, J. A., et al. (2002). "Single-chip surface micromachined integrated gyroscope with 50°/h Allan deviation." *IEEE Journal of Solid-State Circuits*, 37(12), 1860-1866.

### 11.2 Noise Analysis

4. Gabrielson, T. B. (1993). "Mechanical-thermal noise in micromachined acoustic and vibration sensors." *IEEE Transactions on Electron Devices*, 40(5), 903-909.

5. Leland, R. P. (2005). "Mechanical-thermal noise in MEMS gyroscopes." *IEEE Sensors Journal*, 5(3), 493-500.

6. Enz, C. C., & Temes, G. C. (1996). "Circuit techniques for reducing the effects of op-amp imperfections." *Proceedings of the IEEE*, 84(11), 1584-1614.

### 11.3 Random Telegraph Noise

7. Takeuchi, K., et al. (2009). "Understanding random telegraph signal noise in the context of large-scale digital circuits." *IEEE International Electron Devices Meeting*.

8. Tega, N., et al. (2011). "Impact of threshold voltage fluctuation due to random telegraph noise on scaled-down SRAM." *IEEE International Reliability Physics Symposium*.

9. Toledano-Luque, M., et al. (2013). "From mean values to distributions of BTI lifetime of deeply scaled FETs through atomistic understanding of the degradation." *IEEE Transactions on Device and Materials Reliability*, 13(4), 444-450.

### 11.4 Sigma-Delta Modulation

10. Schreier, R., & Temes, G. C. (2004). *Understanding Delta-Sigma Data Converters*. Wiley-IEEE Press.

11. Ortmanns, M., & Gerfers, F. (2006). *Continuous-Time Sigma-Delta A/D Conversion*. Springer.

12. Candy, J. C., & Temes, G. C. (1992). *Oversampling Delta-Sigma Data Converters*. IEEE Press.

### 11.5 Capacitive Sensing

13. Lemkin, M., & Boser, B. E. (1999). "A three-axis micromachined accelerometer with a CMOS position-sense interface and digital offset-trim electronics." *IEEE Journal of Solid-State Circuits*, 34(4), 456-468.

14. Amini, B. V., & Ayazi, F. (2005). "Micro-gravity capacitive silicon-on-insulator accelerometers." *Journal of Micromechanics and Microengineering*, 15(11), 2113.

---

## 12. Appendices

### Appendix A: Complete Noise Source Summary

| Source | Symbol | Location | Type | Typical Value | Section |
|--------|--------|----------|------|---------------|---------|
| **Drive Channel** |
| Drive C/V OTA | $e_{n,opa,drive}$ | Drive C/V input | Thermal + Flicker | 5-10 nV/√Hz | 2.2.1 |
| MEMS Drive Rs | $e_{n,rs,drive}$ | DP/DN traces | Thermal | 1-2 nV/√Hz | 2.2.2 |
| ASIC Drive Rp | $e_{n,rp,drive}$ | DP/DN input | Thermal (current) | 10-30 nV/√Hz @ 35kHz | 2.2.3 |
| Phase Int Rsig | $e_{n,rsig,phase}$ | Phase int input | Thermal | 10-40 nV/√Hz | 2.3.1 |
| Phase Int OTA | $e_{n,opa,phase}$ | Phase int | Thermal + Flicker | 5-10 nV/√Hz | 2.3.2 |
| VCO Phase Noise | $\phi_{vco}$ | VCO output | Colored | -80 to -120 dBc/Hz | 2.3.3 |
| **Sense Channel - References** |
| Bandgap | $e_{n,bg}$ | Reference gen | Thermal + Flicker | 50-100 nV/√Hz | 3.2.1 |
| IDAC | $i_{n,idac}$ | Current DAC | Thermal + Flicker | 0.1-1% FS | 3.2.2 |
| CM LDO | $e_{n,cm}$ | CM voltage | Thermal + Flicker | 100-500 nV/√Hz | 3.2.3 |
| **Sense Channel - MEMS** |
| MEMS Sense Rs | $e_{n,rs,sense}$ | CP/CN traces | Thermal | 1-3 nV/√Hz | 3.3.1 |
| Brownian | $\Omega_{Brw}$ | Proof mass | Mechanical thermal | 0.01-0.1 dps/√Hz | 3.3.2 |
| **Sense Channel - ASIC** |
| Sense C/V OTA | $e_{n,opa,sense}$ | Sense C/V input | Thermal + Flicker + RTN | 5-15 nV/√Hz | 3.4.1 |
| Sense C/V RTN | $e_{n,rtn}$ | Input transistors | RTN | 0.5-3 mV (steps) | 3.4.2 |
| Rate Int Rsig | $e_{n,rsig,rate}$ | Rate int input | Thermal | 10-40 nV/√Hz | 3.5.1 |
| Rate Int OTA | $e_{n,opa,rate}$ | Rate integrator | Thermal + Flicker | 5-10 nV/√Hz | 3.5.2 |
| ASIC Sense Rp | $e_{n,rp,sense}$ | CP/CN input | Thermal (current) | 10-30 nV/√Hz @ 35kHz | 3.6 |
| **CM Coupling Paths** |
| CM over C0 | $e_{n,cm,c0}$ | Capacitor offset | VCM × ΔC0 / C | 5-20 nV/√Hz | 4.2 |
| CM Parasitic | $e_{n,cm,par}$ | Substrate mismatch | VCM × ΔCpar / C | 2-10 nV/√Hz | 4.3 |
| CM Rate Mod | $e_{n,cm,rate}$ | Rate-dependent | VCM × rate | Signal-dependent | 4.4 |
| CM Quad Demod | $e_{n,cm,quad}$ | Quad demodulation | VCM × Qres | 1-5 nV/√Hz | 4.5 |
| QC Coupling | $e_{n,qc}$ | Quad compensation | Drive noise × CQC/C | 2-15 nV/√Hz | 4.6 |
| **Sigma-Delta** |
| SD1 Quantization | $e_{q,sd1}$ | 1st stage quantizer | Shaped quantization | Depends on OSR | 5.2 |
| SD2 Quantization | $e_{q,sd2}$ | 2nd stage quantizer | Shaped quantization | Depends on OSR | 5.3 |
| Phase Noise SD Mod | $e_{n,sd\_mod}$ | Feedback DAC | φ × Qnoise folding | 1-5% increase | 5.2.2 |

### Appendix B: Key Equations Reference

#### Thermal Noise (Johnson-Nyquist)
$$e_n = \sqrt{4 k_B T R} \quad [V/\sqrt{Hz}]$$

#### Brownian Noise (Mechanical)
$$\Omega_{Brw} = \sqrt{\frac{4 k_B T \omega_0}{m_{eff} Q S_{gain}^2 A_{drive}^2}} \quad [dps/\sqrt{Hz}]$$

#### Capacitive Divider Noise Gain
$$G_n = \frac{C_{total}}{C_{feedback}}$$

#### Pseudo-Sinusoidal Noise Coefficient
$$K_{in,noise} = \sqrt{\int_0^T w_{pseudo}^2(t) \, dt} \approx 0.761$$

#### RTN Power Spectral Density (Single Trap)
$$S_{RTN}(f) = \frac{4 A^2 \tau_{avg}}{1 + (2\pi f \tau_{avg})^2} \quad [V^2/Hz]$$

#### Sigma-Delta Noise Transfer Function (MASH 1-1)
$$H_{NTF}(z) = (1 - z^{-1})^2$$

#### Angle Random Walk
$$ARW = \sigma_{white} \times \sqrt{3600} \times 60 \quad [deg/\sqrt{hr}]$$

#### Effective Number of Bits
$$ENOB = \frac{SNR_{dB} - 1.76}{6.02}$$

### Appendix C: Parameter Definitions

**Physical Constants:**
- $k_B = 1.38 \times 10^{-23}$ J/K: Boltzmann constant
- $e = 1.602 \times 10^{-19}$ C: Elementary charge
- $\epsilon_0 = 8.854 \times 10^{-12}$ F/m: Permittivity of free space

**MEMS Parameters:**
- $m_{eff}$: Effective sense mass [kg]
- $Q$: Mechanical quality factor [dimensionless]
- $\omega_0 = 2\pi f_0$: Resonant frequency [rad/s]
- $A_{drive}$: Drive amplitude [m]
- $S_{gain}$: Sensitivity [F/dps] or [m/(dps·m)]

**ASIC Parameters:**
- $C_{CVD}$: Drive C/V feedback capacitor [F]
- $C_{CVS}$: Sense C/V feedback capacitor [F]
- $R_{FB}$: Integrator feedback resistor [Ω]
- $f_s$: Sampling frequency [Hz]

**Noise Descriptors:**
- $e_n$: Voltage noise density [V/√Hz]
- $i_n$: Current noise density [A/√Hz]
- $\Omega_n$: Angular rate noise density [dps/√Hz]
- $\sigma(\tau)$: Allan deviation at averaging time τ [dps or deg/hr]

### Appendix D: MATLAB Utility Functions

#### D.1 Colored Noise Generation

```matlab
function [time, noise] = colouredNoise(duration, fs, freq, psd, seed, ...
                                       offset_enable, offset_value, ...
                                       flicker_enable, verbose)
% COLOUREDNOISE - Generate colored noise with specified PSD
%
% Inputs:
%   duration       - Signal duration [s]
%   fs             - Sampling frequency [Hz]
%   freq           - Frequency points for PSD [Hz]
%   psd            - PSD values at frequency points [unit²/Hz]
%   seed           - Random seed (0 = random)
%   offset_enable  - Add DC offset (0/1)
%   offset_value   - DC offset value [unit]
%   flicker_enable - Enable 1/f noise (0/1)
%   verbose        - Print diagnostics (0/1)
%
% Outputs:
%   time  - Time vector [s]
%   noise - Colored noise signal [unit]
%
% Method:
%   1. Generate white noise in frequency domain
%   2. Shape spectrum according to PSD specification
%   3. IFFT to time domain
%   4. Optionally add 1/f component

    if seed ~= 0
        rng(seed);
    end
    
    N = round(duration * fs);
    time = (0:N-1) / fs;
    
    % Frequency vector
    f = (0:N-1) * fs / N;
    
    % Interpolate PSD to all frequencies
    psd_interp = interp1(freq, psd, f, 'linear', 'extrap');
    psd_interp(psd_interp < 0) = 0;  % No negative PSD
    
    % Generate white noise in frequency domain
    phase = 2*pi*rand(1, N);  % Random phase
    amplitude = sqrt(psd_interp * fs / 2);  % Scale for sampling
    
    % Complex spectrum
    spectrum = amplitude .* exp(1j * phase);
    spectrum(1) = 0;  % Zero DC
    
    % Make Hermitian symmetric for real IFFT
    spectrum(N/2+2:end) = conj(spectrum(N/2:-1:2));
    
    % IFFT to time domain
    noise = real(ifft(spectrum));
    
    % Add 1/f component if requested
    if flicker_enable
        flicker = generate_1f_noise(N, fs);
        noise = noise + flicker * sqrt(mean(psd(freq<100)));
    end
    
    % Add DC offset
    if offset_enable
        noise = noise + offset_value;
    end
    
    % Verbose output
    if verbose
        fprintf('Colored noise generated:\n');
        fprintf('  Duration: %.3f s\n', duration);
        fprintf('  Samples: %d\n', N);
        fprintf('  RMS: %.3e\n', std(noise));
    end
end

function noise_1f = generate_1f_noise(N, fs)
    % Simple 1/f noise generator
    f = (1:N/2) / N * fs;
    amplitude = 1 ./ sqrt(f);
    amplitude(1) = 0;  % Zero DC
    
    phase = 2*pi*rand(1, N/2);
    spectrum = [0, amplitude .* exp(1j * phase)];
    spectrum = [spectrum, conj(spectrum(end:-1:2))];
    
    noise_1f = real(ifft(spectrum));
    noise_1f = noise_1f / std(noise_1f);  % Normalize
end
```

#### D.2 Spectral Analysis Helper

```matlab
function spectrum = plotFunction(signal, varargin)
% PLOTFUNCTION - Comprehensive spectral analysis and plotting
%
% Usage:
%   spectrum = plotFunction(signal, 'fs', 100e3, 'f0', 35e3, ...
%                          'fb', 200, 'plot_fft', true)
%
% Inputs (name-value pairs):
%   signal    - Input signal vector
%   fs        - Sampling frequency [Hz]
%   f0        - Signal frequency (for SNR calculation) [Hz]
%   fb        - Bandwidth for integration [Hz]
%   fsig      - Alternative signal frequency [Hz]
%   stats     - Calculate statistics (true/false)
%   unit      - Signal unit string (e.g., 'dps', 'V')
%   plot_fft  - Generate plot (true/false)
%   new_fig   - Create new figure (true/false)
%
% Output:
%   spectrum - Structure with fields:
%              .psd, .frequency, .rms, .snr, etc.

    % Parse inputs
    p = inputParser;
    addRequired(p, 'signal');
    addParameter(p, 'fs', 100e3);
    addParameter(p, 'f0', 0);
    addParameter(p, 'fb', 200);
    addParameter(p, 'fsig', 0);
    addParameter(p, 'stats', true);
    addParameter(p, 'unit', 'V');
    addParameter(p, 'plot_fft', false);
    addParameter(p, 'new_fig', false);
    parse(p, signal, varargin{:});
    
    % Calculate PSD
    [psd, f] = pwelch(p.Results.signal, [], [], [], p.Results.fs);
    
    % Store results
    spectrum.psd = psd;
    spectrum.frequency = f;
    spectrum.asd = sqrt(psd);
    
    % Calculate statistics
    if p.Results.stats
        % RMS in bandwidth
        idx_band = f <= p.Results.fb;
        spectrum.rms_band = sqrt(sum(psd(idx_band)) * mean(diff(f)));
        
        % Total RMS
        spectrum.rms_total = std(p.Results.signal);
        
        % SNR (if signal frequency specified)
        if p.Results.f0 > 0 || p.Results.fsig > 0
            f_sig = max(p.Results.f0, p.Results.fsig);
            [~, idx_sig] = min(abs(f - f_sig));
            signal_power = psd(idx_sig) * mean(diff(f));
            noise_power = sum(psd([1:idx_sig-2, idx_sig+2:end])) * mean(diff(f));
            spectrum.snr = 10*log10(signal_power / noise_power);
        end
    end
    
    % Plot if requested
    if p.Results.plot_fft
        if p.Results.new_fig
            figure;
        end
        
        loglog(f, sqrt(psd), 'LineWidth', 1.5);
        xlabel('Frequency [Hz]');
        ylabel(sprintf('ASD [%s/√Hz]', p.Results.unit));
        title('Power Spectral Density');
        grid on;
        
        % Add statistics text
        if p.Results.stats
            text(0.02, 0.98, sprintf('RMS: %.3e %s', ...
                 spectrum.rms_total, p.Results.unit), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top');
        end
    end
end
```

### Appendix E: Troubleshooting Guide

#### E.1 Common Issues

**Issue: Total noise much higher than expected**

*Possible Causes:*
1. Incorrect noise source configuration (all enabled when only subset should be)
2. Parameter values incorrect (check resistances, capacitances)
3. Gain calculations wrong (noise gains too high)
4. Missing noise filtering (CM noise not filtered)

*Debug Steps:*
```matlab
% Check individual contributions
report = generate_noise_breakdown(nm);

% Identify dominant source
fprintf('Dominant: %s (%.1f%%)\n', ...
        report.dominant.C.source, report.dominant.C.percentage);

% Disable suspected source and re-run
nm_test = noise_model(imu, sys, init_noise_cfg_enhanced('Disable_CM'));
```

**Issue: Simulation runs very slowly**

*Possible Causes:*
1. Too many samples ($N$ too large)
2. Colored noise generation with many frequency points
3. Multiple tool calls in loop

*Solutions:*
```matlab
% Reduce sample count for testing
obj.ctrl.N = 10000;  % Instead of 1000000

% Reduce frequency points in colored noise
freq_reduced = logspace(log10(1), log10(fs/2), 20);  % Instead of 100 points

% Disable plotting during batch runs
obj.ctrl.plot_spectra_s = 0;
obj.ctrl.plot_spectra_d = 0;
```

**Issue: NaN or Inf in results**

*Possible Causes:*
1. Division by zero (zero capacitance value)
2. Negative value in sqrt() (PSD calculation)
3. Invalid parameter (Inf resistance)

*Debug:*
```matlab
% Check for NaN
if any(isnan(nm.info.C.sense.noise.ratein_r_dps))
    fprintf('NaN detected in rate signal\n');
    
    % Find source
    fields = fieldnames(nm.info.C.sense.noise);
    for i = 1:length(fields)
        if any(isnan(nm.info.C.sense.noise.(fields{i})))
            fprintf('  NaN in: %s\n', fields{i});
        end
    end
end
```

#### E.2 Validation Failures

**Validation vs. SPICE fails (>20% error)**

*Check:*
1. SPICE simulation settings match MATLAB (sampling rate, bandwidth)
2. Temperature matches
3. Component values match exactly
4. Noise models in SPICE are correct (flicker, thermal)

**Validation vs. measurement fails (>30% error)**

*Consider:*
1. Process variations in actual device
2. Temperature during measurement
3. Supply voltage effects
4. Package parasitics not in model
5. PCB noise coupling

#### E.3 Performance Optimization

**For faster simulation:**

```matlab
% Minimal configuration
nm = noise_model(imu, sys, init_noise_cfg_enhanced('Minimal'));

% Parallel processing (if Parallel Computing Toolbox available)
parfor i = 1:N_runs
    nm_array(i) = noise_model(imu, sys, noise_cfg);
end

% Pre-allocate and vectorize
% (Most operations already vectorized, but check custom additions)
```

**For more accurate results:**

```matlab
% Increase sample count
obj.ctrl.N = 1000000;  % 10 seconds @ 100 kHz

% Finer frequency resolution in colored noise
freq_fine = logspace(log10(0.1), log10(fs/2), 200);

% Enable all noise sources
nm = noise_model(imu, sys, init_noise_cfg_enhanced('All'));
```

---

## Document Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2020-2021 | Original Authors | Initial implementation |
| 2.0 | 2024-10 | Ahmed Hussein | Major revision with enhanced documentation |
| 3.0 | 2024-12 | Ahmed Hussein | Comprehensive merge with implementation details, RTN modeling, CM coupling analysis, validation framework |

---

**End of Document**
