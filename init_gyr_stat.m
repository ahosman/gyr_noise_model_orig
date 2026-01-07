function stat = init_gyr_stat(nsim)
    % Initialize gyroscope noise statistics structure for Monte Carlo simulation
    % Returns an array of stat structures, one for each simulation
    % Each stat(ii) has scalar mc values and metadata fields
    
    % Create template structure with explicit field names
    template = struct();
    template.mc = struct('K_f', [], 'i_1Hz', [], 'v_th', []);
    template.spec_min = struct();
    template.spec_max = struct();
    template.name = struct();
    template.mult = struct();
    template.units = struct();
    template.color = struct();
    
    % Replicate template for nsim simulations
    stat = repmat(template, 1, nsim);
    
    % Thermal voltage noise floor v_th [V/sqrt(Hz)]
    % Gaussian distribution, typical range: 5-15 nV/sqrt(Hz)
    v_th_mean = 30.0E-09;    % 300 nV/sqrt(Hz)
    v_th_std = 0.2E-09;      % 2 nV/sqrt(Hz) (sigma)
    
    % Flicker current noise i_1Hz [A/sqrt(Hz)] @ 1Hz
    % Lognormal distribution, typical range: 50-200 aA/sqrt(Hz)
    % Using log-space parameters: ln(120e-18*4) ≈ -32.24
    i_1Hz_mu_ln = log(120.0E-18*4);  % ln(120 aA/sqrt(Hz))
    i_1Hz_sigma_ln = 0.3;        % Log-space sigma
    
    % Populate each simulation
    for ii = 1:nsim
        % Generate v_th sample (Gaussian, clipped to positive)
        v_th_sample = max(v_th_mean + v_th_std * randn(1), 0);
        
        % Generate i_1Hz sample (Lognormal)
        i_1Hz_sample = exp(i_1Hz_mu_ln + i_1Hz_sigma_ln * randn(1));
        
        % Assign to both K_f (legacy) and i_1Hz (new)
        stat(ii).mc.K_f = i_1Hz_sample;
        stat(ii).mc.i_1Hz = i_1Hz_sample;
        stat(ii).mc.v_th = v_th_sample;
        
        % K_f metadata
        stat(ii).spec_min.K_f = NaN;
        stat(ii).spec_max.K_f = NaN;
        stat(ii).name.K_f = 'Flicker Current @1Hz (i_{1Hz})';
        stat(ii).mult.K_f = 1.0e18;  % display in aA/sqrt(Hz)
        stat(ii).units.K_f = 'aA/sqrt(Hz)';
        stat(ii).color.K_f = 'green';
        
        % i_1Hz metadata
        stat(ii).spec_min.i_1Hz = NaN;
        stat(ii).spec_max.i_1Hz = NaN;
        stat(ii).name.i_1Hz = 'Flicker Current @1Hz (i_{1Hz})';
        stat(ii).mult.i_1Hz = 1.0e18;
        stat(ii).units.i_1Hz = 'aA/sqrt(Hz)';
        stat(ii).color.i_1Hz = 'green';
        
        % v_th metadata
        stat(ii).spec_min.v_th = NaN;
        stat(ii).spec_max.v_th = NaN;
        stat(ii).name.v_th = 'Thermal Noise Floor (v_{th})';
        stat(ii).mult.v_th = 1.0e9;  % display in nV/sqrt(Hz)
        stat(ii).units.v_th = 'nV/sqrt(Hz)';
        stat(ii).color.v_th = 'green';
    end
end
