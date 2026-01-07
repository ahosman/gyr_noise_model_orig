function plot_spectra(obj, axis)
% PLOT_SPECTRA - Plot noise spectral density for specified axis
%
% Inputs:
%   obj   - noise_model object
%   axis  - 'C', 'S', 'Z', or 'all' (default: 'all')
%
% Usage:
%   nm = noise_model(imu, sys, 'All');
%   nm.plot_spectra('all');    % All three axes
%   nm.plot_spectra('C');      % Single axis

if nargin < 2
    axis = 'all';
end

if strcmpi(axis, 'all')
    figure('Name', 'Noise Spectral Density - All Axes', 'NumberTitle', 'off');
    
    for i = 1:3
        ax_name = {'C', 'S', 'Z'};
        subplot(3,1,i);
        
        signal = obj.info.(ax_name{i}).sense.noise.ratein_r_dps;
        [psd, f] = pwelch(signal, [], [], [], obj.imu.as.gyr.config.inf.fs);
        
        loglog(f, sqrt(psd)*1000, 'LineWidth', 1.5);
        xlabel('Frequency [Hz]');
        ylabel('ASD [mdps/√Hz]');
        title(sprintf('%s-Axis Noise Spectral Density', ax_name{i}));
        grid on;
        xlim([1 obj.ctrl.fb]);
        
        % Add RMS annotation
        rms_val = obj.get_total_noise(ax_name{i});
        text(0.05, 0.95, sprintf('RMS: %.3f mdps', rms_val*1000), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'w', 'FontSize', 10);
    end
else
    if ~ismember(upper(axis), {'C', 'S', 'Z'})
        error('Invalid axis: %s (use C, S, Z, or all)', axis);
    end
    
    figure('Name', sprintf('Noise Spectral Density - %s Axis', axis), 'NumberTitle', 'off');
    
    signal = obj.info.(axis).sense.noise.ratein_r_dps;
    [psd, f] = pwelch(signal, [], [], [], obj.imu.as.gyr.config.inf.fs);
    
    loglog(f, sqrt(psd)*1000, 'LineWidth', 2);
    xlabel('Frequency [Hz]');
    ylabel('ASD [mdps/√Hz]');
    title(sprintf('%s-Axis Noise Spectral Density', axis));
    grid on;
    xlim([1 obj.ctrl.fb]);
    
    % Add text with RMS value and statistics
    rms_val = obj.get_total_noise(axis);
    text(0.05, 0.95, sprintf('RMS: %.3f mdps', rms_val*1000), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w', 'FontSize', 11);
end
