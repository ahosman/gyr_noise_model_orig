function [as_sns_cv, gyr_trm_cap_prog_cv_C, gyr_trm_cap_prog_cv_S, gyr_trm_cap_prog_cv_Z] = init_sense_cv(varargin)
phys=phys_pkg.phys;

asic_struct = varargin{1};
asic_config = varargin{2};
if length(varargin) > 2
    gyro_struct = varargin{3};
end
as_sns_cv = struct;
if length(varargin) > 2
    as_sns_cv.cv_out_target    = 1.2;
    as_sns_cv.cv_gain_target   = as_sns_cv.cv_out_target/asic_config.inf.FSdps; % [V/dps]
    
    sense_fe_cv_cap_cv = asic_struct.sense.fe.cv.cap_cv;
    if asic_config.otp.ana.gyr_trm_cvs_fs_2x_en == 1
        sense_fe_cv_cap_cv(3) = [];
        sense_fe_cv_cap_cv(1) = [];
    else
        sense_fe_cv_cap_cv(8) = [];
        sense_fe_cv_cap_cv(2) = [];
    end
    
    %% estimation capacitor values for sense loop
    gyr_cap_cv_s_v = zeros(1,length(2^(length(sense_fe_cv_cap_cv)-1)));
    for ii = 1:(2^(length(sense_fe_cv_cap_cv)-1)-1)
        gyr_cap_cv_s_v(ii) =  sum(sense_fe_cv_cap_cv.*dec2binvec(ii-1, 5, 0));
    end
    
    as_sns_cv.surface_charge       = 0.0;%gyro_struct.mc_results.process.surface_voltage;
    
    %% trimming of sense detection
    detection_target_v = (gyr_cap_cv_s_v*as_sns_cv.cv_gain_target) / (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge);
    [res_error_C, index_C] = min(abs(detection_target_v./gyro_struct.stat.mc.Sdd_C_Fdps - 1));
    [res_error_S, index_S] = min(abs(detection_target_v./gyro_struct.stat.mc.Sdd_S_Fdps - 1));
    [res_error_Z, index_Z] = min(abs(detection_target_v./gyro_struct.stat.mc.Sdd_Z_Fdps - 1));
    gyr_trm_cap_prog_cv_C = index_C - 1;
    gyr_trm_cap_prog_cv_S = index_S - 1;
    gyr_trm_cap_prog_cv_Z = index_Z - 1;
    
    as_sns_cv.C.gyr_cap_cv_s = sum(sense_fe_cv_cap_cv.*dec2binvec(gyr_trm_cap_prog_cv_C, 5, 0));
    as_sns_cv.S.gyr_cap_cv_s = sum(sense_fe_cv_cap_cv.*dec2binvec(gyr_trm_cap_prog_cv_S, 5, 0));
    as_sns_cv.Z.gyr_cap_cv_s = sum(sense_fe_cv_cap_cv.*dec2binvec(gyr_trm_cap_prog_cv_Z, 5, 0));
    
    %% C/V ouput at FSdps (should be 800mV, single ended peak to peak or zero-peak fully differential for FS=6000dps)
    if asic_config.otp.ana.gyr_trm_cvs_fs_2x_en == 1
        as_sns_cv.C.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_C_Fdps*asic_config.inf.FSdps ...
            /(2.0*as_sns_cv.C.gyr_cap_cv_s);
        as_sns_cv.S.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_S_Fdps*asic_config.inf.FSdps ...
            /(2.0*as_sns_cv.S.gyr_cap_cv_s);
        as_sns_cv.Z.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_Z_Fdps*asic_config.inf.FSdps ...
            /(2.0*as_sns_cv.Z.gyr_cap_cv_s);
    else
        as_sns_cv.C.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_C_Fdps*asic_config.inf.FSdps ...
            /(as_sns_cv.C.gyr_cap_cv_s);
        as_sns_cv.S.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_S_Fdps*asic_config.inf.FSdps ...
            /(as_sns_cv.S.gyr_cap_cv_s);
        as_sns_cv.Z.cv_out = ...
            (asic_config.inf.ref.cm.vout - as_sns_cv.surface_charge)*gyro_struct.stat.mc.Sdd_Z_Fdps*asic_config.inf.FSdps ...
            /(as_sns_cv.Z.gyr_cap_cv_s);
    end
    
    %% Gain C/V [V/dps]   -> should be 2e-4 V/dps
    as_sns_cv.C.cv_gain = as_sns_cv.C.cv_out / asic_config.inf.FSdps;
    as_sns_cv.S.cv_gain = as_sns_cv.S.cv_out / asic_config.inf.FSdps;
    as_sns_cv.Z.cv_gain = as_sns_cv.Z.cv_out / asic_config.inf.FSdps;
    
    %% Gain C/V from Iin to Vout [V/A]
    as_sns_cv.C.cv_gain_Iin_Vout = 1 / (2*pi*asic_config.inf.drive.vco.gyro_drv_fres*as_sns_cv.C.gyr_cap_cv_s);
    as_sns_cv.S.cv_gain_Iin_Vout = 1 / (2*pi*asic_config.inf.drive.vco.gyro_drv_fres*as_sns_cv.S.gyr_cap_cv_s);
    as_sns_cv.Z.cv_gain_Iin_Vout = 1 / (2*pi*asic_config.inf.drive.vco.gyro_drv_fres*as_sns_cv.Z.gyr_cap_cv_s);
    
    as_sns_cv.C.gyr_dgain = round( res_error_C / (0.25/2^9));
    as_sns_cv.S.gyr_dgain = round( res_error_S / (0.25/2^9));
    as_sns_cv.Z.gyr_dgain = round( res_error_Z / (0.25/2^9));
    
    as_sns_cv.C.cin_pair  = 300.0E-15;
    as_sns_cv.S.cin_pair  = 300.0E-15;
    as_sns_cv.Z.cin_pair  = 300.0E-15;
    
    %% Noise
    switch asic_config.def.corner
        case 'T'
            % C/V amplifier noise [V/rt(Hz)]
            if strcmp(asic_config.def.mode,'HPM')
                as_sns_cv.noise.freq = [1.000 1.080 1.166 1.259 1.359 1.468 1.585 1.711 1.848 1.995 2.154 2.326 2.512 2.712 2.929 3.162 3.415 3.687 3.981 4.299 4.642 5.012 5.412 5.843 6.310 6.813 7.356 7.943 8.577 9.261 10.00 ...
                    10.80 11.66 12.59 13.59 14.68 15.85 17.11 18.48 19.95 21.54 23.26 25.12 27.12 29.29 31.62 34.15 36.87 39.81 42.99 46.42 50.12 54.12 58.43 63.10 68.13 73.56 79.43 85.77 92.61 100.0 ...
                    108.0 116.6 125.9 135.9 146.8 158.5 171.1 184.8 199.5 215.4 232.6 251.2 271.2 292.9 316.2 341.5 368.7 398.1 429.9 464.2 501.2 541.2 584.3 631.0 681.3 735.6 794.3 857.7 926.1 1.000E3 ...
                    1.080E3 1.166E3 1.259E3 1.359E3 1.468E3 1.585E3 1.711E3 1.848E3 1.995E3 2.154E3 2.326E3 2.512E3 2.712E3 2.929E3 3.162E3 3.415E3 3.687E3 3.981E3 4.299E3 4.642E3 5.012E3 5.412E3 5.843E3 6.310E3 6.813E3 7.356E3 7.943E3 8.577E3 9.261E3 10.00E3 ...
                    10.80E3 11.66E3 12.59E3 13.59E3 14.68E3 15.85E3 17.11E3 18.48E3 19.95E3 21.54E3 23.26E3 25.00E3 25.12E3 27.12E3 29.29E3 31.62E3 34.15E3 35.00E3 36.87E3 39.81E3 42.99E3 46.42E3 50.12E3 54.12E3 58.43E3 63.10E3 68.13E3 73.56E3 79.43E3 85.77E3 92.61E3 100.0E3 ...
                    108.0E3 116.6E3 125.9E3 135.9E3 146.8E3 158.5E3 171.1E3 184.8E3 199.5E3 215.4E3 232.6E3 251.2E3 271.2E3 292.9E3 316.2E3 341.5E3 368.7E3 398.1E3 429.9E3 464.2E3 501.2E3 541.2E3 584.3E3 631.0E3 681.3E3 735.6E3 794.3E3 857.7E3 926.1E3 1.000E6 ...
                    1.080E6 1.166E6 1.259E6 1.359E6 1.468E6 1.585E6 1.711E6 1.848E6 1.995E6 2.154E6 2.326E6 2.512E6 2.712E6 2.929E6 3.162E6 3.415E6 3.687E6 3.981E6 4.299E6 4.642E6 5.012E6 5.412E6 5.843E6 6.310E6 6.813E6 7.356E6 7.943E6 8.577E6 9.261E6 10.00E6 ...
                    10.80E6 11.66E6 12.59E6 13.59E6 14.68E6 15.85E6 17.11E6 18.48E6 19.95E6 21.54E6 23.26E6 25.12E6 27.12E6 29.29E6 31.62E6 34.15E6 36.87E6 39.81E6 42.99E6 46.42E6 50.12E6 54.12E6 58.43E6 63.10E6 68.13E6 73.56E6 79.43E6 85.77E6 92.61E6 100.0E6 ...
                    108.0E6 116.6E6 125.9E6 135.9E6 146.8E6 158.5E6 171.1E6 184.8E6 199.5E6 215.4E6 232.6E6 251.2E6 271.2E6 292.9E6 316.2E6 341.5E6 368.7E6 398.1E6 429.9E6 464.2E6 501.2E6 541.2E6 584.3E6 631.0E6 681.3E6 735.6E6 794.3E6 857.7E6 926.1E6 1.000E9];
                % as_sns_cv.noise.vnd = [8880.0 8235.0 7638.0 7085.0 6573.0 6098.0 5658.0 5251.0 4874.0 4525.0 ...
                %     4201.0 3902.0 3624.0 3367.0 3128.0 2908.0 2703.0 2514.0 2338.0 2176.0 ...
                %     2025.0 1885.0 1756.0 1636.0 1525.0 1422.0 1326.0 1238.0 1156.0 1080.0 ...
                %     1009.0 943.7 883.0 826.6 774.3 725.7 680.6 638.7 599.7 563.5 529.9 498.6 ...
                %     469.4 442.3 417.0 393.5 371.5 351.0 331.9 314.0 297.3 281.7 267.1 253.4 ...
                %     240.6 228.6 217.3 206.7 196.7 187.3 178.5 170.2 162.3 154.9 147.9 141.3 ...
                %     135.1 129.1 123.5 118.2 113.2 108.4 103.9 99.55 95.45 91.54 87.83 84.29 ...
                %     80.92 77.71 74.65 71.72 68.94 66.28 63.74 61.31 59.00 56.79 54.67 52.65 ...
                %     50.72 48.88 47.12 45.43 43.82 42.28 40.81 39.40 38.06 36.77 35.55 34.37 ...
                %     33.25 32.19 31.17 30.19 29.26 28.38 27.54 26.73 25.97 25.24 24.54 23.88 ...
                %     23.26 22.66 22.09 21.56 21.05 20.56 20.11 19.67 19.26 18.88 18.51 18.16 ...
                %     17.84 17.53 17.24 16.96 16.71 16.46 16.25 16.24 16.02 15.82 15.63 15.46 ...
                %     15.40 15.29 15.13 14.99 14.85 14.72 14.60 14.49 14.39 14.29 14.20 14.11 ...
                %     14.03 13.96 13.89 13.83 13.77 13.71 13.66 13.61 13.57 13.53 13.49 13.45 ...
                %     13.42 13.39 13.37 13.34 13.32 13.30 13.28 13.27 13.26 13.25 13.25 13.25 ...
                %     13.25 13.26 13.27 13.29 13.32 13.36 13.40 13.46 13.53 13.61 13.71 13.82 ...
                %     13.95 14.10 14.26 14.42 14.60 14.76 14.92 15.07 15.19 15.29 15.37 15.42 ...
                %     15.46 15.47 15.48 15.47 15.45 15.43 15.41 15.39 15.37 15.35 15.34 15.33 ...
                %     15.32 15.32 15.33 15.33 15.35 15.37 15.40 15.43 15.47 15.52 15.58 15.64 ...
                %     15.73 15.82 15.93 16.07 16.23 16.43 16.67 16.97 17.33 17.77 18.31 18.95 ...
                %     19.73 20.65 21.76 23.07 24.62 26.45 28.58 31.07 33.95 37.28 41.10 45.45 ...
                %     50.39 55.93 62.08 68.80 75.95 83.32 90.54 97.15 102.7 106.7 109.0 109.6 ...
                %     108.8 106.9 104.2 101.1 97.78 94.34 90.88 87.44 84.05 80.69 77.37 74.07 70.79 67.52 64.27]*1.0E-9;
                v_th = asic_struct.stat.mc.v_th;% V/√Hz
                K_f = asic_struct.stat.mc.K_f;% V/√Hz @1Hz
                
                f_ref = 1;
                gamma = 1.0; % coeff
                % as_sns_cv.noise.vnd = v_th*sqrt(phys.T/300)*ones(1,length(as_sns_cv.noise.freq));
                as_sns_cv.noise.vnd = sqrt(v_th^2*ones(1,length(as_sns_cv.noise.freq))+K_f^2*f_ref./as_sns_cv.noise.freq.^(2.0*gamma));
                as_sns_cv.noise.Vn_opa = 23.8e-9*sqrt(phys.T/300);
            end
            if strcmp(asic_config.def.mode,'NORM')
                as_sns_cv.noise.Vn_opa = 32e-9*sqrt(phys.T/300);   % C/V amplifier noise in low power mode [V/rt(Hz)]
            end
            if strcmp(asic_config.def.mode,'gm0')
                as_sns_cv.noise.Vn_opa = 18e-9*sqrt(phys.T/300);
            end
            if strcmp(asic_config.def.mode, 'gm1')
                as_sns_cv.noise.Vn_opa = 23e-9*sqrt(phys.T/300);
            end
        case 'S'
            as_sns_cv.noise.Vn_opa = 16e-9*sqrt(phys.T/300); %17.1e-9*sqrt(phys.T/300);   % C/V amplifier noise [V/rt(Hz)]
            if asic_config.def.lp_en == 1
                as_sns_cv.noise.Vn_opa = as_sns_cv.noise.Vn_opa * sqrt((3/4) / (1/4));   % C/V amplifier noise in low power mode [V/rt(Hz)]
            end
        case 'F'
            as_sns_cv.noise.Vn_opa = 14e-9*sqrt(phys.T/300); %15.3e-9*sqrt(phys.T/300);   % C/V amplifier noise [V/rt(Hz)]
            if asic_config.def.lp_en == 1
                as_sns_cv.noise.Vn_opa = as_sns_cv.noise.Vn_opa * sqrt((3/4) / (1/4));   % C/V amplifier noise in low power mode [V/rt(Hz)]
            end
    end
end
end
