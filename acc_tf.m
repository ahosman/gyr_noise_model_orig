% MATLAB Program for Accelerometer Transfer Function and ASIC Requirements

% This script calculates the transfer function of a 3-axis MEMS capacitive-type accelerometer,
% including trimming and digital calibration. It also defines the requirements for various parts
% of the ASIC, including optimal register sizes and resolutions, based on expected variations
% and limits of parameters from the MEMS sensor.

% Author: Osman, Ahmed Hussein (BST/EIC6)
% Date: 04.10.2024

function acc_tf(obj)

    % Accelerometer Specifications
    full_scale_range = 8;        % Using ±8g range
    digital_resolution = 24;     % bits
    ADC_full_scale = 2;          % Assuming ±1 V full-scale range
    ADC_resolution = digital_resolution;  % bits
    
    % Physical Constants
    g = 9.80665;  % acceleration due to gravity in m/s^2
    
    %% 2. Define Expected Variations and Limits from MEMS
    
    % Temperature Range
    T_min = -40;  % in °C
    T_max = 85;   % in °C
    delta_T = T_max - 25;  % From reference temperature to max temp
    
    % Zero-g Offset Variations
    O_T = 0.8e-3;  % Zero-g offset temp coefficient in g/K
    delta_O_temp = O_T * delta_T;  % Offset variation due to temperature
    
    offset_soldering_max = 50e-3;  % in g (±50 mg)
    
    strain_max = 50;  % in µstrain (Assumed maximum strain)
    O_strain_coeff = 0.02e-3;  % Zero-g offset change with respect to bending (g/µstrain)
    delta_O_bending = O_strain_coeff * strain_max;  % Offset variation due to bending
    
    delta_O_total_g = abs(delta_O_temp) + offset_soldering_max + delta_O_bending;
    
    % Sensitivity Variations
    sensitivity_error_life = 0.75 / 100;  % ±0.75%
    alpha_S = 0.01 / 100;  % Sensitivity temp coefficient (%/K)
    delta_S_temp = alpha_S * delta_T;  % Sensitivity variation due to temperature
    
    beta_S = 0.08 / 100;  % Sensitivity drift over bending (%/µstrain)
    delta_S_bending = beta_S * strain_max;  % Sensitivity variation due to bending
    
    sensitivity_nonlinearity_temp = 0.1 / 100;  % ±0.1%
    sensitivity_nonlinearity_8g = 1 / 100;  % Max 1%
    sensitivity_error_8g = 1 / 100;  % Max 1%
    
    sensitivity_variations = [
        sensitivity_error_life;
        delta_S_temp;
        delta_S_bending;
        sensitivity_nonlinearity_temp;
        sensitivity_nonlinearity_8g;
        sensitivity_error_8g
    ];
    
    total_sensitivity_error_percent = sqrt(sum(sensitivity_variations .^ 2)) * 100;  % in %
    
    %% 3. Calculate Transfer Function Including Trimming and Calibration
    
    FS_output = 2^(digital_resolution - 1) - 1;  % for a bipolar output
    S_nominal = FS_output / full_scale_range;    % counts per g
    
    %% 4. Determine Optimal Register Sizes and Resolutions
    
    % Offset Correction Register Size
    delta_O_total_counts = delta_O_total_g * S_nominal;
    offset_corr_bits = ceil(log2(delta_O_total_counts + 1)) + 1;  % +1 for sign bit
    
    fprintf('Total Zero-g Offset Variation: %.3f mg\n', delta_O_total_g * 1e3);
    fprintf('Offset Correction in Counts: %.0f counts\n', delta_O_total_counts);
    fprintf('Offset Correction Register Size: %d bits\n', offset_corr_bits);
    
    % Gain Correction Register Size
    total_sensitivity_error = total_sensitivity_error_percent;  % in %
    gain_correction_resolution_percent = 0.001;  % 0.001% resolution
    gain_adjustment_range_percent = total_sensitivity_error * 2;  % ± total error
    gain_steps = gain_adjustment_range_percent / gain_correction_resolution_percent;
    gain_corr_bits = ceil(log2(gain_steps + 1)) + 1;  % +1 for sign bit
    
    fprintf('Total Sensitivity Error: ±%.2f%%\n', total_sensitivity_error);
    fprintf('Gain Correction Steps Needed: %.0f\n', gain_steps);
    fprintf('Gain Correction Register Size: %d bits\n', gain_corr_bits);
    
    %% 5. Output the Results
    
    fprintf('\n--- ASIC Requirements ---\n');
    fprintf('Offset Correction Register Size: %d bits\n', offset_corr_bits);
    fprintf('Gain Correction Register Size: %d bits\n', gain_corr_bits);
    
    fprintf('\nTransfer Function:\n');
    fprintf('  E_corrected = [E_raw - O_corr(T, sigma)] * G_corr(T, sigma)\n');

end
