function report = generate_report(obj, options)
if nargin < 2
    options = struct();
end

% Set defaults
if ~isfield(options, 'signal_band')
    options.signal_band = 200;
end
if ~isfield(options, 'plot')
    options.plot = true;
end
if ~isfield(options, 'verbose')
    options.verbose = true;
end
if ~isfield(options, 'save_report')
    options.save_report = false;
end
if ~isfield(options, 'filename')
    options.filename = 'noise_report.txt';
end

% Generate report structure
report = struct();
report.config = obj.config;
report.options = options;

% Noise summary for each axis
for axis = {'C', 'S', 'Z'}
    ax = axis{1};
    report.(ax).rms_total = obj.get_total_noise(ax);
    report.(ax).rms_mdps = report.(ax).rms_total * 1000;
    
    % Spectral analysis
    signal = obj.info.(ax).sense.noise.ratein_r_dps;
    [psd, f] = pwelch(signal, [], [], [], obj.imu.as.gyr.config.inf.fs);
    
    % RMS in signal band
    idx_band = f <= options.signal_band;
    report.(ax).rms_in_band = sqrt(trapz(f(idx_band), psd(idx_band)));
    
    % Store frequency and PSD for detailed analysis
    report.(ax).frequency = f;
    report.(ax).psd = psd;
    
    % Calculate peak frequency
    [~, idx_peak] = max(psd);
    report.(ax).peak_freq = f(idx_peak);
    report.(ax).peak_value = sqrt(psd(idx_peak));
end

% Print report to console
if options.verbose
    obj.print_noise_report(report);
end

% Save report to file
if options.save_report
    obj.save_report_file(report, options.filename);
end
