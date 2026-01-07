function stat = init_gyro_stat(nsim, txt_tbl_5000)
% Drive Mass
stat.drive.drive_mass = struct;
stat.drive.drive_mass.mean     = 1.2859e-08; %(gyro.mc_results.drive.modal_drive_displacement_for_1_m)^2;
stat.drive.drive_mass.std      = 1.007e-9;
stat.drive.drive_mass.min      = 18.6354e-9;
stat.drive.drive_mass.max      = 21.1874e-9;
stat.drive.drive_mass.spec_min = 18.6354e-9;
stat.drive.drive_mass.spec_max = 21.1874e-9;

stat.sense.sense_mass_C = struct;
stat.sense.sense_mass_C.mean = 4.2011e-09;
stat.sense.sense_mass_C.std  = 1.007e-9;
stat.sense.sense_mass_C.min  = 0.8*stat.sense.sense_mass_C.mean;
stat.sense.sense_mass_C.max  = 1.2*stat.sense.sense_mass_C.mean;
stat.sense.sense_mass_C.spec_min  = 0.8*stat.sense.sense_mass_C.mean;
stat.sense.sense_mass_C.spec_max  = 1.2*stat.sense.sense_mass_C.mean;

stat.sense.sense_mass_S = struct;
stat.sense.sense_mass_S.mean = 6.4648e-09;
stat.sense.sense_mass_S.std  = 1.007e-9;
stat.sense.sense_mass_S.min  = 0.8*stat.sense.sense_mass_S.mean;
stat.sense.sense_mass_S.max  = 1.2*stat.sense.sense_mass_S.mean;
stat.sense.sense_mass_S.spec_min  = 0.8*stat.sense.sense_mass_S.mean;
stat.sense.sense_mass_S.spec_max  = 1.2*stat.sense.sense_mass_S.mean;

stat.sense.sense_mass_Z = struct;
stat.sense.sense_mass_Z.mean = 5.9328e-09;
stat.sense.sense_mass_Z.std  = 1.007e-9;
stat.sense.sense_mass_Z.min  = 0.8*stat.sense.sense_mass_Z.mean;
stat.sense.sense_mass_Z.max  = 1.2*stat.sense.sense_mass_Z.mean;
stat.sense.sense_mass_Z.spec_min  = 0.8*stat.sense.sense_mass_Z.mean;
stat.sense.sense_mass_Z.spec_max  = 1.2*stat.sense.sense_mass_Z.mean;

% Use default field values (Excel dependency retired)
stat.drive.quality_factor  = create_default_field(); % Drive Quality Factor
stat.drive.drive_amplitude = create_default_field(); % Drive Amplitude in um
stat.drive.sensitivity.aa  = create_default_field(); % AA Actuation Sensitivity in F/m
stat.drive.sensitivity.ai  = create_default_field(); % AI Actuation Sensitivity in F/m
stat.drive.sensitivity.da  = create_default_field(); % DA Detection Sensitivity in F/m
stat.drive.sensitivity.di  = create_default_field(); % DI Detection Sensitivity in F/m
stat.drive.capacitance.aa  = create_default_field(); % AA Actuation capacitance
stat.drive.capacitance.ai  = create_default_field();
stat.drive.capacitance.da  = create_default_field();
stat.drive.capacitance.di  = create_default_field();

% Drive Frequency in Hz
stat.drive.tuned_frequency = create_default_field();

% Gyro Sense Quality Factor
stat.sense.C.quality_factor = create_default_field();
stat.sense.S.quality_factor = create_default_field();
stat.sense.Z.quality_factor = create_default_field();

% Gyro Sense Tuned Frequency
stat.sense.C.tuned_frequency = create_default_field();
stat.sense.S.tuned_frequency = create_default_field();
stat.sense.Z.tuned_frequency = create_default_field();

% Gyro Sense Quadrature
stat.sense.C.quadrature = create_default_field();
stat.sense.S.quadrature = create_default_field();
stat.sense.Z.quadrature = create_default_field();

stat.sense.C.capacitance.cp = create_default_field();
stat.sense.C.capacitance.cn = create_default_field();

stat.sense.S.capacitance.cp = create_default_field();
stat.sense.S.capacitance.cn = create_default_field();

stat.sense.Z.capacitance.cp = create_default_field();
stat.sense.Z.capacitance.cn = create_default_field();

% Gyro Sense Quadrature Drift
stat.sense.C.quadrature_drift = struct;
stat.sense.C.quadrature_drift.mean     = 0.0;
stat.sense.C.quadrature_drift.std      = 1;
stat.sense.C.quadrature_drift.min      = -1500;
stat.sense.C.quadrature_drift.max      = +1500;
stat.sense.C.quadrature_drift.spec_min = -1500;
stat.sense.C.quadrature_drift.spec_max = +1500;

stat.sense.S.quadrature_drift = struct;
stat.sense.S.quadrature_drift.mean     = 0.0;
stat.sense.S.quadrature_drift.std      = 1;
stat.sense.S.quadrature_drift.min      = -1500;
stat.sense.S.quadrature_drift.max      = +1500;
stat.sense.S.quadrature_drift.spec_min = -1500;
stat.sense.S.quadrature_drift.spec_max = +1500;

stat.sense.Z.quadrature_drift = struct;
stat.sense.Z.quadrature_drift.mean     = 0.0;
stat.sense.Z.quadrature_drift.std      = 1;
stat.sense.Z.quadrature_drift.min      = -1500;
stat.sense.Z.quadrature_drift.max      = +1500;
stat.sense.Z.quadrature_drift.spec_min = -1500;
stat.sense.Z.quadrature_drift.spec_max = +1500;

% Gyro Sense Sensitivity
stat.sense.C.sensitivity = create_default_field();
stat.sense.S.sensitivity = create_default_field();
stat.sense.Z.sensitivity = create_default_field();

% Gyro tuned frequency split
stat.sense.C.tuned_split = create_default_field();
stat.sense.S.tuned_split = create_default_field();
stat.sense.Z.tuned_split = create_default_field();

stat.process.surface_voltage  = 0.0;%gyro.process.surface_voltage;

% Create MC samples from txt_tbl_5000
stat_tmp.mc = create_sample(stat, txt_tbl_5000, nsim);

if true % Always execute MC mode logic
    % Save the base template structure and replicate it for all iterations
    stat_base = stat;
    stat = repmat(stat_base, 1, nsim);
    
    for ii = 1:nsim
        stat(ii).mc.xd                 = stat_tmp.mc.xd(ii);
        stat(ii).spec_min.xd           = stat(ii).drive.drive_amplitude.min;
        stat(ii).spec_max.xd           = stat(ii).drive.drive_amplitude.max;
        stat(ii).name.xd               = "Drive Amplitude (xd)";
        stat(ii).mult.xd               = 1.0E+06;
        stat(ii).units.xd              = 'um';
        stat(ii).color.xd              = 'blue';
        
        stat(ii).mc.md                 = stat_tmp.mc.md(ii);
        stat(ii).spec_min.md           = stat(ii).drive.drive_mass.spec_min;
        stat(ii).spec_max.md           = stat(ii).drive.drive_mass.spec_max;
        stat(ii).name.md               = "Drive Mass (md)";
        stat(ii).mult.md               = 1.0;
        stat(ii).units.md              = 'kg';
        stat(ii).color.md              = 'blue';
        
        stat(ii).mc.Qd                 = stat_tmp.mc.Qd(ii);
        stat(ii).spec_min.Qd           = stat(ii).drive.quality_factor.min;
        stat(ii).spec_max.Qd           = stat(ii).drive.quality_factor.max;
        stat(ii).name.Qd               = "Drive Quality Factor (Qd)";
        stat(ii).mult.Qd               = 1.0;
        stat(ii).units.Qd              = ' ';
        stat(ii).color.Qd              = 'blue';
        
        stat(ii).mc.fd                 = stat_tmp.mc.fd(ii);
        stat(ii).spec_min.fd           = stat(ii).drive.tuned_frequency.spec_min;
        stat(ii).spec_max.fd           = stat(ii).drive.tuned_frequency.spec_max;
        stat(ii).name.fd               = "Drive Resonance Frequency (fd)";
        stat(ii).mult.fd               = 1.0E-03;
        stat(ii).units.fd              = 'kHz';
        stat(ii).color.fd              = 'blue';
        
        stat(ii).mc.Saa                = stat_tmp.mc.Saa(ii);
        stat(ii).spec_min.Saa          = stat(ii).drive.sensitivity.aa.min;
        stat(ii).spec_max.Saa          = stat(ii).drive.sensitivity.aa.max;
        stat(ii).name.Saa              = "Drive Actuation Sensitivity (Saa)";
        stat(ii).mult.Saa              = 1.0E+9;
        stat(ii).units.Saa             = 'fF/um';
        stat(ii).color.Saa             = 'blue';
        
        stat(ii).mc.Sai                = stat_tmp.mc.Sai(ii);
        stat(ii).spec_min.Sai          = stat(ii).drive.sensitivity.ai.min;
        stat(ii).spec_max.Sai          = stat(ii).drive.sensitivity.ai.max;
        stat(ii).name.Sai              = "Drive Actuation Sensitivity (Sai)";
        stat(ii).mult.Sai              = 1.0E+9;
        stat(ii).units.Sai             = 'fF/um';
        stat(ii).color.Sai             = 'blue';
        
        stat(ii).mc.Sda                = stat_tmp.mc.Sda(ii);
        stat(ii).spec_min.Sda          = stat(ii).drive.sensitivity.da.min;
        stat(ii).spec_max.Sda          = stat(ii).drive.sensitivity.da.max;
        stat(ii).name.Sda              = "Drive Detection Sensitivity (Sda)";
        stat(ii).mult.Sda              = 1.0E+9;
        stat(ii).units.Sda             = 'fF/um';
        stat(ii).color.Sda             = 'blue';
        
        stat(ii).mc.Sdi                = stat_tmp.mc.Sdi(ii);
        stat(ii).spec_min.Sdi          = stat(ii).drive.sensitivity.di.min;
        stat(ii).spec_max.Sdi          = stat(ii).drive.sensitivity.di.max;
        stat(ii).name.Sdi              = "Drive Detection Sensitivity (Sdi)";
        stat(ii).mult.Sdi              = 1.0E+9;
        stat(ii).units.Sdi             = 'fF/um';
        stat(ii).color.Sdi             = 'blue';
        
        stat(ii).mc.Cap                = stat_tmp.mc.Cap(ii);
        stat(ii).spec_min.Cap          = stat(ii).drive.capacitance.aa.min;
        stat(ii).spec_max.Cap          = stat(ii).drive.capacitance.aa.max;
        stat(ii).name.Cap              = "Drive Actuation Working Capacitance (Cap)";
        stat(ii).mult.Cap              = 1.0E+15;
        stat(ii).units.Cap             = 'fF';
        stat(ii).color.Cap             = 'blue';
        
        stat(ii).mc.Can                = stat_tmp.mc.Can(ii);
        stat(ii).spec_min.Can          = stat(ii).drive.capacitance.ai.min;
        stat(ii).spec_max.Can          = stat(ii).drive.capacitance.ai.max;
        stat(ii).name.Can              = "Drive Actuation Working Capacitance (Can)";
        stat(ii).mult.Can              = 1.0E+15;
        stat(ii).units.Can             = 'fF';
        stat(ii).color.Can             = 'blue';
        
        stat(ii).mc.Cdp                = stat_tmp.mc.Cdp(ii);
        stat(ii).spec_min.Cdp          = stat(ii).drive.capacitance.da.min;
        stat(ii).spec_max.Cdp          = stat(ii).drive.capacitance.da.max;
        stat(ii).name.Cdp              = "Drive Detection Working Capacitance (Cdp)";
        stat(ii).mult.Cdp              = 1.0E+15;
        stat(ii).units.Cdp             = 'fF';
        stat(ii).color.Cdp             = 'blue';
        
        stat(ii).mc.Cdn                = stat_tmp.mc.Cdn(ii);
        stat(ii).spec_min.Cdn          = stat(ii).drive.capacitance.di.min;
        stat(ii).spec_max.Cdn          = stat(ii).drive.capacitance.di.max;
        stat(ii).name.Cdn              = "Drive Detection Working Capacitance (Cdn)";
        stat(ii).mult.Cdn              = 1.0E+15;
        stat(ii).units.Cdn             = 'fF';
        stat(ii).color.Cdn             = 'blue';
        
        stat(ii).mc.Csp_C              = stat_tmp.mc.Csp_C(ii);
        stat(ii).spec_min.Csp_C        = stat(ii).sense.C.capacitance.cp.min;
        stat(ii).spec_max.Csp_C        = stat(ii).sense.C.capacitance.cp.max;
        stat(ii).name.Csp_C            = "Sense Working Capacitance - C (Csp\_C)";
        stat(ii).mult.Csp_C            = 1.0E+15;
        stat(ii).units.Csp_C           = 'fF';
        stat(ii).color.Csp_C           = 'blue';
        
        stat(ii).mc.Csn_C              = stat_tmp.mc.Csn_C(ii);
        stat(ii).spec_min.Csn_C        = stat(ii).sense.C.capacitance.cn.min;
        stat(ii).spec_max.Csn_C        = stat(ii).sense.C.capacitance.cn.max;
        stat(ii).name.Csn_C            = "Sense Working Capacitance - C (Csn\_C)";
        stat(ii).mult.Csn_C            = 1.0E+15;
        stat(ii).units.Csn_C           = 'fF';
        stat(ii).color.Csn_C           = 'blue';
        
        stat(ii).mc.Csp_S              = stat_tmp.mc.Csp_S(ii);
        stat(ii).spec_min.Csp_S        = stat(ii).sense.S.capacitance.cp.min;
        stat(ii).spec_max.Csp_S        = stat(ii).sense.S.capacitance.cp.max;
        stat(ii).name.Csp_S            = "Sense Working Capacitance - S (Csp\_S)";
        stat(ii).mult.Csp_S            = 1.0E+15;
        stat(ii).units.Csp_S           = 'fF';
        stat(ii).color.Csp_S           = 'blue';
        
        stat(ii).mc.Csn_S              = stat_tmp.mc.Csn_S(ii);
        stat(ii).spec_min.Csn_S        = stat(ii).sense.S.capacitance.cn.min;
        stat(ii).spec_max.Csn_S        = stat(ii).sense.S.capacitance.cn.max;
        stat(ii).name.Csn_S            = "Sense Working Capacitance - S (Csn\_S)";
        stat(ii).mult.Csn_S            = 1.0E+15;
        stat(ii).units.Csn_S           = 'fF';
        stat(ii).color.Csn_S           = 'blue';
        
        stat(ii).mc.Csp_Z              = stat_tmp.mc.Csp_Z(ii);
        stat(ii).spec_min.Csp_Z        = stat(ii).sense.Z.capacitance.cp.min;
        stat(ii).spec_max.Csp_Z        = stat(ii).sense.Z.capacitance.cp.max;
        stat(ii).name.Csp_Z            = "Sense Working Capacitance - Z (Csp\_Z)";
        stat(ii).mult.Csp_Z            = 1.0E+15;
        stat(ii).units.Csp_Z           = 'fF';
        stat(ii).color.Csp_Z           = 'blue';
        
        stat(ii).mc.Csn_Z              = stat_tmp.mc.Csn_Z(ii);
        stat(ii).spec_min.Csn_Z        = stat(ii).sense.Z.capacitance.cn.min;
        stat(ii).spec_max.Csn_Z        = stat(ii).sense.Z.capacitance.cn.max;
        stat(ii).name.Csn_Z            = "Sense Working Capacitance - Z (Csn\_Z)";
        stat(ii).mult.Csn_Z            = 1.0E+15;
        stat(ii).units.Csn_Z           = 'fF';
        stat(ii).color.Csn_Z           = 'blue';
        
        stat(ii).mc.ms_C               = stat_tmp.mc.ms_C(ii);
        stat(ii).spec_min.ms_C         = stat(ii).sense.sense_mass_C.spec_min;
        stat(ii).spec_max.ms_C         = stat(ii).sense.sense_mass_C.spec_max;
        stat(ii).name.ms_C             = "Sense Mass C (ms\_C)";
        stat(ii).mult.ms_C             = 1.0;
        stat(ii).units.ms_C            = 'kg';
        stat(ii).color.ms_C            = 'blue';
        
        stat(ii).mc.ms_S               = stat_tmp.mc.ms_S(ii);
        stat(ii).spec_min.ms_S         = stat(ii).sense.sense_mass_S.spec_min;
        stat(ii).spec_max.ms_S         = stat(ii).sense.sense_mass_S.spec_max;
        stat(ii).name.ms_S             = "Sense Mass S (ms\_S)";
        stat(ii).mult.ms_S             = 1.0;
        stat(ii).units.ms_S            = 'kg';
        stat(ii).color.ms_S            = 'blue';
        
        stat(ii).mc.ms_Z               = stat_tmp.mc.ms_Z(ii);
        stat(ii).spec_min.ms_Z         = stat(ii).sense.sense_mass_Z.spec_min;
        stat(ii).spec_max.ms_Z         = stat(ii).sense.sense_mass_Z.spec_max;
        stat(ii).name.ms_Z             = "Sense Mass Z (ms\_Z)";
        stat(ii).mult.ms_Z             = 1.0;
        stat(ii).units.ms_Z            = 'kg';
        stat(ii).color.ms_Z            = 'blue';
        
        stat(ii).mc.Qs_C               = stat_tmp.mc.Qs_C(ii);
        stat(ii).spec_min.Qs_C         = stat(ii).sense.C.quality_factor.spec_min;
        stat(ii).spec_max.Qs_C         = stat(ii).sense.C.quality_factor.spec_max;
        stat(ii).name.Qs_C             = "Sense Quality Factor -C (Qs\_C)";
        stat(ii).mult.Qs_C             = 1.0;
        stat(ii).units.Qs_C            = ' ';
        stat(ii).color.Qs_C            = 'blue';
        
        stat(ii).mc.Qs_S               = stat_tmp.mc.Qs_S(ii);
        stat(ii).spec_min.Qs_S         = stat(ii).sense.S.quality_factor.spec_min;
        stat(ii).spec_max.Qs_S         = stat(ii).sense.S.quality_factor.spec_max;
        stat(ii).name.Qs_S             = "Sense Quality Factor -S (Qs\_S)";
        stat(ii).mult.Qs_S             = 1.0;
        stat(ii).units.Qs_S            = ' ';
        stat(ii).color.Qs_S            = 'blue';
        
        stat(ii).mc.Qs_Z               = stat_tmp.mc.Qs_Z(ii);
        stat(ii).spec_min.Qs_Z         = stat(ii).sense.Z.quality_factor.spec_min;
        stat(ii).spec_max.Qs_Z         = stat(ii).sense.Z.quality_factor.spec_max;
        stat(ii).name.Qs_Z             = "Sense Quality Factor -Z (Qs\_Z)";
        stat(ii).mult.Qs_Z             = 1.0;
        stat(ii).units.Qs_Z            = ' ';
        stat(ii).color.Qs_Z            = 'blue';
        
        stat(ii).mc.fs_C               = stat_tmp.mc.fs_C(ii);
        stat(ii).spec_min.fs_C         = stat(ii).sense.C.tuned_frequency.min;
        stat(ii).spec_max.fs_C         = stat(ii).sense.C.tuned_frequency.max;
        stat(ii).name.fs_C             = "Sense Resonance Frequency -C (fs\_C)";
        stat(ii).mult.fs_C             = 1.0E-03;
        stat(ii).units.fs_C            = 'kHz';
        stat(ii).color.fs_C            = 'blue';
        
        stat(ii).mc.fs_S               = stat_tmp.mc.fs_S(ii);
        stat(ii).spec_min.fs_S         = stat(ii).sense.S.tuned_frequency.min;
        stat(ii).spec_max.fs_S         = stat(ii).sense.S.tuned_frequency.max;
        stat(ii).name.fs_S             = "Sense Resonance Frequency -S (fs\_S)";
        stat(ii).mult.fs_S             = 1.0E-03;
        stat(ii).units.fs_S            = 'kHz';
        stat(ii).color.fs_S            = 'blue';
        
        stat(ii).mc.fs_Z               = stat_tmp.mc.fs_Z(ii);
        stat(ii).spec_min.fs_Z         = stat(ii).sense.Z.tuned_frequency.min;
        stat(ii).spec_max.fs_Z         = stat(ii).sense.Z.tuned_frequency.max;
        stat(ii).name.fs_Z             = "Sense Resonance Frequency -Z (fs\_Z)";
        stat(ii).mult.fs_Z             = 1.0E-03;
        stat(ii).units.fs_Z            = 'kHz';
        stat(ii).color.fs_Z            = 'blue';
        
        stat(ii).mc.fsplit_C           = stat_tmp.mc.fsplit_C(ii);
        stat(ii).spec_min.fsplit_C     = stat(ii).sense.C.tuned_split.spec_min;
        stat(ii).spec_max.fsplit_C     = stat(ii).sense.C.tuned_split.spec_max;
        stat(ii).name.fsplit_C         = "Sense Frequency Split -C (\Delta f\_C)";
        stat(ii).mult.fsplit_C         = 1.0E-03;
        stat(ii).units.fsplit_C        = 'kHz';
        stat(ii).color.fsplit_C        = 'blue';
        
        stat(ii).mc.fsplit_S           = stat_tmp.mc.fsplit_S(ii);
        stat(ii).spec_min.fsplit_S     = stat(ii).sense.S.tuned_split.spec_min;
        stat(ii).spec_max.fsplit_S     = stat(ii).sense.S.tuned_split.spec_max;
        stat(ii).name.fsplit_S         = "Sense Frequency Split -S (\Delta f\_S)";
        stat(ii).mult.fsplit_S         = 1.0E-03;
        stat(ii).units.fsplit_S        = 'kHz';
        stat(ii).color.fsplit_S        = 'blue';
        
        stat(ii).mc.fsplit_Z           = stat_tmp.mc.fsplit_Z(ii);
        stat(ii).spec_min.fsplit_Z     = stat(ii).sense.Z.tuned_split.spec_min;
        stat(ii).spec_max.fsplit_Z     = stat(ii).sense.Z.tuned_split.spec_max;
        stat(ii).name.fsplit_Z         = "Sense Frequency Split -Z (\Delta f\_Z)";
        stat(ii).mult.fsplit_Z         = 1.0E-03;
        stat(ii).units.fsplit_Z        = 'kHz';
        stat(ii).color.fsplit_Z        = 'blue';
        
        stat(ii).mc.quad_C             = stat_tmp.mc.quad_C(ii);
        stat(ii).spec_min.quad_C       = stat(ii).sense.C.quadrature.min;
        stat(ii).spec_max.quad_C       = stat(ii).sense.C.quadrature.max;
        stat(ii).name.quad_C           = "Quadrature -C (quad\_C)";
        stat(ii).mult.quad_C           = 1.0;
        stat(ii).units.quad_C          = 'dps';
        stat(ii).color.quad_C          = 'blue';
        
        stat(ii).mc.quad_S             = stat_tmp.mc.quad_S(ii);
        stat(ii).spec_min.quad_S       = stat(ii).sense.S.quadrature.min;
        stat(ii).spec_max.quad_S       = stat(ii).sense.S.quadrature.max;
        stat(ii).name.quad_S           = "Quadrature -S (quad\_S)";
        stat(ii).mult.quad_S           = 1.0;
        stat(ii).units.quad_S          = 'dps';
        stat(ii).color.quad_S          = 'blue';
        
        stat(ii).mc.quad_Z             = stat_tmp.mc.quad_Z(ii);
        stat(ii).spec_min.quad_Z       = stat(ii).sense.Z.quadrature.min;
        stat(ii).spec_max.quad_Z       = stat(ii).sense.Z.quadrature.max;
        stat(ii).name.quad_Z           = "Quadrature -Z (quad\_Z)";
        stat(ii).mult.quad_Z           = 1.0;
        stat(ii).units.quad_Z          = 'dps';
        stat(ii).color.quad_Z          = 'blue';
        
        stat(ii).mc.quad_drift_C       = stat_tmp.mc.quad_drift_C(ii);
        stat(ii).spec_min.quad_drift_C = stat(ii).sense.C.quadrature_drift.spec_min;
        stat(ii).spec_max.quad_drift_C = stat(ii).sense.C.quadrature_drift.spec_max;
        stat(ii).name.quad_drift_C     = "Quadrature Drift -C (quad\_drift\_C)";
        stat(ii).mult.quad_drift_C     = 1.0;
        stat(ii).units.quad_drift_C    = 'dps';
        stat(ii).color.quad_drift_C    = 'blue';
        
        stat(ii).mc.quad_drift_S       = stat_tmp.mc.quad_drift_S(ii);
        stat(ii).spec_min.quad_drift_S = stat(ii).sense.S.quadrature_drift.spec_min;
        stat(ii).spec_max.quad_drift_S = stat(ii).sense.S.quadrature_drift.spec_max;
        stat(ii).name.quad_drift_S     = "Quadrature Drift -S (quad\_drift\_S)";
        stat(ii).mult.quad_drift_S     = 1.0;
        stat(ii).units.quad_drift_S    = 'dps';
        stat(ii).color.quad_drift_S    = 'blue';
        
        stat(ii).mc.quad_drift_Z       = stat_tmp.mc.quad_drift_Z(ii);
        stat(ii).spec_min.quad_drift_Z = stat(ii).sense.Z.quadrature_drift.spec_min;
        stat(ii).spec_max.quad_drift_Z = stat(ii).sense.Z.quadrature_drift.spec_max;
        stat(ii).name.quad_drift_Z     = "Quadrature Drift -Z (quad\_drift\_Z)";
        stat(ii).mult.quad_drift_Z     = 1.0;
        stat(ii).units.quad_drift_Z    = 'dps';
        stat(ii).color.quad_drift_Z    = 'blue';
        
        stat(ii).mc.Sdd_C_Fdps              = stat_tmp.mc.Sdd_C_Fdps(ii);
        stat(ii).spec_min.Sdd_C_Fdps        = stat(ii).sense.C.sensitivity.spec_min;
        stat(ii).spec_max.Sdd_C_Fdps        = stat(ii).sense.C.sensitivity.spec_max;
        stat(ii).name.Sdd_C_Fdps            = "Electrical Sensitivitiy -C (Sdd\_C)";
        stat(ii).mult.Sdd_C_Fdps            = 1.0E+18;
        stat(ii).units.Sdd_C_Fdps           = 'aF/dps';
        stat(ii).color.Sdd_C_Fdps           = 'blue';
        
        stat(ii).mc.Sdd_S_Fdps              = stat_tmp.mc.Sdd_S_Fdps(ii);
        stat(ii).spec_min.Sdd_S_Fdps        = stat(ii).sense.S.sensitivity.spec_min;
        stat(ii).spec_max.Sdd_S_Fdps        = stat(ii).sense.S.sensitivity.spec_max;
        stat(ii).name.Sdd_S_Fdps            = "Electrical Sensitivitiy -S (Sdd\_S)";
        stat(ii).mult.Sdd_S_Fdps            = 1.0E+18;
        stat(ii).units.Sdd_S_Fdps           = 'aF/dps';
        stat(ii).color.Sdd_S_Fdps           = 'blue';
        
        stat(ii).mc.Sdd_Z_Fdps              = stat_tmp.mc.Sdd_Z_Fdps(ii);
        stat(ii).spec_min.Sdd_Z_Fdps        = stat(ii).sense.Z.sensitivity.spec_min;
        stat(ii).spec_max.Sdd_Z_Fdps        = stat(ii).sense.Z.sensitivity.spec_max;
        stat(ii).name.Sdd_Z_Fdps            = "Electrical Sensitivitiy -Z (Sdd\_Z)";
        stat(ii).mult.Sdd_Z_Fdps            = 1.0E+18;
        stat(ii).units.Sdd_Z_Fdps           = 'aF/dps';
        stat(ii).color.Sdd_Z_Fdps           = 'blue';
    end
end
end

function sample = create_sample(stat, txt_tbl_5000, nsim)
    % Always read from txt_tbl_5000 (mode parameter retired)
    opts = detectImportOptions(txt_tbl_5000);
    opts.VariableNamingRule = 'preserve';
    tbl_5000            = readtable(txt_tbl_5000, opts);
    sample.xd           = tbl_5000.FP_drive_amplitude_drive;
    sample.md           = tbl_5000.FP_misc_equiv_mass_drive;
    sample.Qd           = tbl_5000.FP_drive_drive_quality_factor;
    sample.fd           = tbl_5000.FP_working_point_tuned_frequency_with_drive_drive;
    
    sample.Sda          = tbl_5000.FP_drive_sensitivity_d1;
    sample.Sdi          = tbl_5000.FP_drive_sensitivity_d2;
    
    sample.Saa          = tbl_5000.FP_drive_sensitivity_a1;
    sample.Sai          = tbl_5000.FP_drive_sensitivity_a2;
    
    sample.Cap          = tbl_5000.FP_drive_capacitance_a1;
    sample.Can          = tbl_5000.FP_drive_capacitance_a2;
    
    sample.Cdp          = tbl_5000.FP_drive_capacitance_d1;
    sample.Cdn          = tbl_5000.FP_drive_capacitance_d2;
    
    sample.ms_C         = tbl_5000.FP_misc_equiv_mass_C;
    sample.ms_S         = tbl_5000.FP_misc_equiv_mass_S;
    sample.ms_Z         = tbl_5000.FP_misc_equiv_mass_Z;
    
    sample.Qs_C         = tbl_5000.FP_detection_quality_factor_C;
    sample.Qs_S         = tbl_5000.FP_detection_quality_factor_S;
    sample.Qs_Z         = tbl_5000.FP_detection_quality_factor_Z;
    
    sample.fs_C         = tbl_5000.FP_working_point_tuned_frequency_with_drive_C;
    sample.fs_S         = tbl_5000.FP_working_point_tuned_frequency_with_drive_S;
    sample.fs_Z         = tbl_5000.FP_working_point_tuned_frequency_with_drive_Z;
    
    sample.fsplit_C     = tbl_5000.FP_working_point_tuned_split_with_drive_C;
    sample.fsplit_S     = tbl_5000.FP_working_point_tuned_split_with_drive_S;
    sample.fsplit_Z     = tbl_5000.FP_working_point_tuned_split_with_drive_Z;
    
    sample.quad_C       = tbl_5000.FP_detection_quadrature_C;
    sample.quad_S       = tbl_5000.FP_detection_quadrature_S;
    sample.quad_Z       = tbl_5000.FP_detection_quadrature_Z;
    
    sample.quad_drift_C = stat.sense.C.quadrature_drift.mean + stat.sense.C.quadrature_drift.std*randn(1,nsim);
    sample.quad_drift_S = stat.sense.S.quadrature_drift.mean + stat.sense.S.quadrature_drift.std*randn(1,nsim);
    sample.quad_drift_Z = stat.sense.Z.quadrature_drift.mean + stat.sense.Z.quadrature_drift.std*randn(1,nsim);
    
    sample.Sdd_C_Fdps        = tbl_5000.FP_detection_corrected_electrical_sensitivity_C;
    sample.Sdd_S_Fdps        = tbl_5000.FP_detection_corrected_electrical_sensitivity_S;
    sample.Sdd_Z_Fdps        = tbl_5000.FP_detection_corrected_electrical_sensitivity_Z;
    
    sample.Csp_C        = tbl_5000.FP_detection_working_point_capacitance_C_cp;
    sample.Csn_C        = tbl_5000.FP_detection_working_point_capacitance_C_cn;
    
    sample.Csp_S        = tbl_5000.FP_detection_working_point_capacitance_S_cp;
    sample.Csn_S        = tbl_5000.FP_detection_working_point_capacitance_S_cn;
    
    sample.Csp_Z        = tbl_5000.FP_detection_working_point_capacitance_Z_cp;
    sample.Csn_Z        = tbl_5000.FP_detection_working_point_capacitance_Z_cn;
    
    
    sample.brw_noise_C  = tbl_5000.FP_misc_brownian_noise_C;
    sample.brw_noise_S  = tbl_5000.FP_misc_brownian_noise_S;
    sample.brw_noise_Z  = tbl_5000.FP_misc_brownian_noise_Z;
end

function field = create_default_field()
    % Create default field structure with nominal values
    % (Excel dependency retired - using defaults)
    field = struct;
    field.mean     = 0.0;
    field.std      = 0.0;
    field.min      = 0.0;
    field.max      = 0.0;
    field.spec_min = 0.0;
    field.spec_max = 0.0;
end
