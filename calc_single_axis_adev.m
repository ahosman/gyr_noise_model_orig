function adev_result = calc_single_axis_adev(obj, signal)
% CALC_SINGLE_AXIS_ADEV - Calculate Allan deviation for single axis
%
% Inputs:
%   obj    - noise_model object
%   signal - Noise signal time series [dps]
%
% Output:
%   adev_result - Structure with:
%               .tau  - Time constant vector [seconds]
%               .adev - Allan deviation values [dps]
%
% Private method - called by allandev()
%
% Usage: (internal use only)
%   adev = obj.calc_single_axis_adev(signal);

% Sample sizes to analyze
adev_result.tau = logspace(-6, 0, 100);  % Time constant [seconds]
adev_result.adev = zeros(size(adev_result.tau));

% Calculate Allan variance for each tau
dt = 1 / obj.imu.as.gyr.config.inf.fs;

for i = 1:length(adev_result.tau)
    tau = adev_result.tau(i);
    N_samples = max(1, round(tau / dt));
    
    if N_samples < length(signal)
        % Integrate signal over tau windows
        N_windows = floor(length(signal) / N_samples);
        windows = reshape(signal(1:N_windows*N_samples), N_samples, N_windows);
        integrated = sum(windows, 1) * dt;
        
        % Calculate variance of integrated samples (Allan variance)
        % ADEV = sqrt(0.5 * mean((theta_n+1 - theta_n)^2)) / (2*tau)
        if N_windows > 1
            differences = diff(integrated);
            adev_result.adev(i) = sqrt(mean(differences.^2)) / (2*sqrt(2)*tau);
        else
            adev_result.adev(i) = NaN;
        end
    else
        adev_result.adev(i) = NaN;
    end
end
