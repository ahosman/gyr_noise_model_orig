function [adev_C, adev_S, adev_Z] = allandev(obj)
% ALLANDEV - Calculate Allan deviation for all three axes
%
% Outputs:
%   adev_C - Allan deviation structure for C-axis
%   adev_S - Allan deviation structure for S-axis
%   adev_Z - Allan deviation structure for Z-axis
%
% Each output structure contains:
%   .tau  - Time constant vector [seconds]
%   .adev - Allan deviation values [dps]
%
% Note: This can be time-consuming for large datasets
%
% Usage:
%   nm = noise_model(imu, sys, 'All');
%   nm.ctrl.calc_adev = 1;
%   [adev_c, adev_s, adev_z] = nm.allandev();

fprintf('Calculating Allan deviation...\n');

% Calculate for each axis
adev_C = obj.calc_single_axis_adev(obj.info.C.sense.noise.ratein_r_dps);
adev_S = obj.calc_single_axis_adev(obj.info.S.sense.noise.ratein_r_dps);
adev_Z = obj.calc_single_axis_adev(obj.info.Z.sense.noise.ratein_r_dps);

fprintf('  ✓ C-axis: %d time constants\n', length(adev_C.tau));
fprintf('  ✓ S-axis: %d time constants\n', length(adev_S.tau));
fprintf('  ✓ Z-axis: %d time constants\n', length(adev_Z.tau));
