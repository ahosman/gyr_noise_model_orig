function w = pseudo_sin(t, N)
% PSEUDO_SIN - Generate pseudo-sinusoidal waveform
%
% This function generates a stepped approximation to a sinusoidal waveform
% used for demodulation in MEMS gyroscope readout circuits.
%
% Inputs:
%   t - Time indices (0, 1, 2, ..., length-1)
%   N - Period in samples (typically 192)
%
% Output:
%   w - Pseudo-sinusoidal waveform with values:
%       +1.0, +1/2.4, 0, -1/2.4, -1.0
%
% Waveform Properties:
%   Period: N samples (192 typical at 100 kHz → ~520 Hz)
%   Levels: 5 discrete levels
%   RMS value: ~0.761 (noise coefficient K_in_noise)
%   THD: ~15-20% compared to pure sine
%
% Phase Assignments (for N=192):
%   [4-24]:    +1/2.4 (rising edge)
%   [24-72]:   +1.0   (peak)
%   [72-94]:   +1/2.4 (falling edge)
%   [94-98]:    0     (zero crossing)
%   [98-120]:  -1/2.4 (falling edge)
%   [120-168]: -1.0   (trough)
%   [168-190]: -1/2.4 (rising edge)
%   [0-4, 190-192]: 0 (zero crossing)
%
% References: noise_model.md Section 6
%
% Example:
%   w = noise_model.pseudo_sin(0:191, 192);
%   plot(w);
%
% See also: noise_model

% Vectorized implementation for efficiency
frac = mod(t, N) ./ N;  % Normalized phase [0, 1)
w = zeros(size(t));

% Define middle level (1/2.4 = 0.4167)
mid_level = 1/2.4;

% Positive half-cycle
w(frac >= 4/192 & frac < 24/192) = mid_level;    % Rising edge
w(frac >= 24/192 & frac < 72/192) = 1.0;         % Peak
w(frac >= 72/192 & frac < 94/192) = mid_level;   % Falling edge

% Zero crossing at 94-98/192 (already initialized to 0)

% Negative half-cycle
w(frac >= 98/192 & frac < 120/192) = -mid_level;  % Falling edge
w(frac >= 120/192 & frac < 168/192) = -1.0;       % Trough
w(frac >= 168/192 & frac < 190/192) = -mid_level; % Rising edge

% Zero crossing at 0-4/192 and 190-192 (already initialized to 0)
