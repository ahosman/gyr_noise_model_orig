function print_noise_report(obj, report)
% PRINT_NOISE_REPORT - Print formatted noise report to console
%
% Inputs:
%   obj    - noise_model object
%   report - Report structure from generate_report()
%
% Private method - called by generate_report()
%
% Usage: (internal use only)
%   obj.print_noise_report(report);

fprintf('\n');
fprintf('=======================================================\n');
fprintf('                   NOISE ANALYSIS REPORT\n');
fprintf('=======================================================\n');
fprintf('Configuration: %s\n', report.config.mode);
fprintf('Timestamp: %s\n', report.config.timestamp);
fprintf('Sample rate: %.0f kHz\n', report.config.fs/1000);
fprintf('Duration: %.2f seconds\n', report.config.duration);
fprintf('Total samples: %d\n', report.config.N);
fprintf('\n');

fprintf('Total RMS Noise by Axis:\n');
fprintf('-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf('  %s-Axis: %10.4f dps  (%8.3f mdps)\n', ax, ...
            report.(ax).rms_total, report.(ax).rms_mdps);
end
fprintf('\n');

fprintf('RMS Noise in Signal Band (0-%d Hz):\n', report.options.signal_band);
fprintf('-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf('  %s-Axis: %10.6f dps  (%8.4f mdps)\n', ax, ...
            report.(ax).rms_in_band, report.(ax).rms_in_band*1000);
end
fprintf('\n');

fprintf('Spectral Characteristics:\n');
fprintf('-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf('  %s-Axis:\n', ax);
    fprintf('    Peak frequency: %8.1f Hz\n', report.(ax).peak_freq);
    fprintf('    Peak ASD:       %8.4f mdps/√Hz\n', report.(ax).peak_value*1000);
end
fprintf('\n');

fprintf('=======================================================\n\n');
