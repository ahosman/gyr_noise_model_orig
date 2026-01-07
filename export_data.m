function export_data(obj, filename)
% EXPORT_DATA - Export noise data to MAT file
%
% Inputs:
%   obj      - noise_model object
%   filename - Output filename (default: 'noise_model_data.mat')
%
% Usage:
%   nm = noise_model(imu, sys, 'All');
%   nm.export_data('my_results.mat');

if nargin < 2
    filename = 'noise_model_data.mat';
end

% Create export structure
export_data = struct();
export_data.config = obj.config;
export_data.info = obj.info;
export_data.spectrum = obj.spectrum;

if obj.ctrl.calc_adev
    export_data.adev = obj.adev;
end

% Save to file
save(filename, 'export_data', '-v7.3');
fprintf('Data exported to: %s\n', filename);
fprintf('  Config:   %s\n', obj.config.mode);
fprintf('  Duration: %.2f s (%d samples)\n', obj.config.duration, obj.config.N);
fprintf('  File size: %.1f MB\n', dir(filename).bytes / 1e6);
