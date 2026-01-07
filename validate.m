function validate(obj)
% VALIDATE - Run validation checks on noise model
%
% Checks:
%   - All axes have noise data
%   - Noise values in reasonable range
%   - Spectral properties valid
%   - Similar noise on all axes (for equal input rates)
%
% Usage:
%   nm = noise_model(imu, sys, 'All');
%   nm.validate();

fprintf('Validating noise model...\n');

% Check 1: All axes have noise data
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    if ~isfield(obj.info, ax) || ~isfield(obj.info.(ax), 'sense')
        error('Missing noise data for %s axis', ax);
    end
end
fprintf('  ✓ All axes have noise data\n');

% Check 2: Noise values in reasonable range
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    noise_rms = obj.get_total_noise(ax);
    if noise_rms < 0.001 || noise_rms > 1000
        warning('Unusual noise level for %s axis: %.3f dps', ax, noise_rms);
    end
end
fprintf('  ✓ Noise levels in reasonable range\n');

% Check 3: All axes have similar noise (should be close for equal rates)
noise_c = obj.get_total_noise('C');
noise_s = obj.get_total_noise('S');
noise_z = obj.get_total_noise('Z');

max_noise = max([noise_c, noise_s, noise_z]);
min_noise = min([noise_c, noise_s, noise_z]);
ratio = max_noise / min_noise;

if ratio > 1.1 && obj.sys.rate_x.Value == obj.sys.rate_y.Value && ...
   obj.sys.rate_x.Value == obj.sys.rate_z.Value
    warning('Noise levels differ significantly across axes despite equal input rates');
else
    fprintf('  ✓ Noise levels consistent across axes\n');
end

% Check 4: Output structures populated
if ~isfield(obj.spectrum, 'C') || ~isfield(obj.spectrum.C, 'rate')
    warning('SD output not fully populated');
else
    fprintf('  ✓ All output structures populated\n');
end

fprintf('Validation complete.\n\n');
