function save_report_file(obj, report, filename)
% SAVE_REPORT_FILE - Save noise report to text file
%
% Inputs:
%   obj      - noise_model object
%   report   - Report structure from generate_report()
%   filename - Output filename
%
% Private method - called by generate_report()
%
% Usage: (internal use only)
%   obj.save_report_file(report, 'noise_report.txt');

fid = fopen(filename, 'w');

if fid < 0
    error('Could not open file for writing: %s', filename);
end

fprintf(fid, '=======================================================\n');
fprintf(fid, '                   NOISE ANALYSIS REPORT\n');
fprintf(fid, '=======================================================\n');
fprintf(fid, 'Configuration: %s\n', report.config.mode);
fprintf(fid, 'Timestamp: %s\n', report.config.timestamp);
fprintf(fid, 'Sample rate: %.0f kHz\n', report.config.fs/1000);
fprintf(fid, 'Duration: %.2f seconds\n', report.config.duration);
fprintf(fid, 'Total samples: %d\n', report.config.N);
fprintf(fid, '\n');

fprintf(fid, 'Total RMS Noise by Axis:\n');
fprintf(fid, '-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf(fid, '  %s-Axis: %10.4f dps  (%8.3f mdps)\n', ax, ...
            report.(ax).rms_total, report.(ax).rms_mdps);
end
fprintf(fid, '\n');

fprintf(fid, 'RMS Noise in Signal Band (0-%d Hz):\n', report.options.signal_band);
fprintf(fid, '-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf(fid, '  %s-Axis: %10.6f dps  (%8.4f mdps)\n', ax, ...
            report.(ax).rms_in_band, report.(ax).rms_in_band*1000);
end
fprintf(fid, '\n');

fprintf(fid, 'Spectral Characteristics:\n');
fprintf(fid, '-------------------------------------------------------\n');
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    fprintf(fid, '  %s-Axis:\n', ax);
    fprintf(fid, '    Peak frequency: %8.1f Hz\n', report.(ax).peak_freq);
    fprintf(fid, '    Peak ASD:       %8.4f mdps/√Hz\n', report.(ax).peak_value*1000);
end
fprintf(fid, '\n');

fprintf(fid, '=======================================================\n');

fclose(fid);
fprintf('Report saved to: %s\n', filename);
