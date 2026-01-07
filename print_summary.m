function print_summary(obj)
% PRINT_SUMMARY - Print quick summary of noise results
%
% Usage:
%   nm = noise_model(imu, sys, 'All');
%   nm.print_summary();

fprintf('\n--- Total Noise Summary ---\n');
fprintf('Configuration: %s\n', obj.config.mode);
fprintf('Duration: %.2f s (%d samples @ %.0f kHz)\n', ...
        obj.config.duration, obj.config.N, obj.config.fs/1000);
fprintf('\n');

fprintf('Total RMS Noise by Axis:\n');
fprintf('-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    rms_dps = obj.get_total_noise(ax);
    rms_mdps = rms_dps * 1000;
    fprintf('  %s-axis: %10.4f dps  (%8.3f mdps)\n', ax, rms_dps, rms_mdps);
end

fprintf('\n');
