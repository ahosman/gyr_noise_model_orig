function [gyr_drv, as_drv_cv, cap_prog_cv_d] = init_drive_cv(varargin)
    phys=phys_pkg.phys;
    
    asic_struct = varargin{1};
    asic_config = varargin{2};
    if length(varargin) > 2
        gyro_struct = varargin{3};
    end
    drive_fe_cv_cap_cv = asic_struct.drive.fe.cv.cap_cv;
    if asic_config.otp.ana.gyr_trm_drv_4um_en == 1
        gyr_drv.drive_ampl_target = 4e-6;
        drive_fe_cv_cap_cv(7) = [];
        drive_fe_cv_cap_cv(1) = [];
    else
        gyr_drv.drive_ampl_target = 8e-6;
        drive_fe_cv_cap_cv(8) = [];
        drive_fe_cv_cap_cv(2) = [];
    end
    gyr_cap_cv_d_v = zeros(1,length(2^(length(drive_fe_cv_cap_cv)-1)));

    for ii = 1:(2^(length(drive_fe_cv_cap_cv)-1))
        gyr_cap_cv_d_v(ii) =  sum(drive_fe_cv_cap_cv.*dec2binvec(ii-1, 5, 0));
    end

    if length(varargin) > 2
        if asic_config.otp.ana.gyr_trm_drv_4um_en == 1
            as_drv_cv.detection_target     = ...
                gyro_struct.stat.mc.Sda*gyro_struct.stat.mc.xd/2.0;
        else
            as_drv_cv.detection_target     = ...
                gyro_struct.stat.mc.Sda*gyro_struct.stat.mc.xd;
        end
        as_drv_cv.drive_ampl = as_drv_cv.detection_target/gyro_struct.stat.mc.Sda;
    else
        % drive detection target single ended signal
        as_drv_cv.gyr_cap_cv_d = sum(drive_fe_cv_cap_cv.*dec2binvec(asic_config.otp.ana.gyr_trm_cap_prog_cv_d, 5, 0));
        as_drv_cv.detection_target = (as_drv_cv.gyr_cap_cv_d*asic_config.inf.drive.agc.out) / ...
            (2*asic_config.inf.ref.cm.vout);
    end

    detection_target_v = (gyr_cap_cv_d_v*asic_config.inf.drive.agc.out) / ...
        (2*asic_config.inf.ref.cm.vout);

    % find new gyr_trm_cap_prog_cv_d for optimal amplitude trimming
    [~, index] = min(abs(detection_target_v./as_drv_cv.detection_target - 1));
    cap_prog_cv_d = index - 1;
    as_drv_cv.gyr_cap_cv_d = sum(drive_fe_cv_cap_cv.*dec2binvec(cap_prog_cv_d, 5, 0));

    if length(varargin) > 2
        as_drv_cv.vout = gyro_struct.stat.mc.Sda*gyro_struct.stat.mc.xd*...
            asic_config.inf.ref.cm.vout/as_drv_cv.gyr_cap_cv_d;
    end
    
    %% Noise Contribution
    
    switch asic_config.def.corner
        case 'T'
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 2)
                as_drv_cv.noise.Vn_opa = 36.3E-09*sqrt(phys.T/300);
            end
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 3)
                as_drv_cv.noise.Vn_opa = 50.8E-09*sqrt(phys.T/300);
            end
        case 'S'
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 2)
                as_drv_cv.noise.Vn_opa = 27.6e-9*sqrt(phys.T/300);
            end
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 3)
                as_drv_cv.noise.Vn_opa = 35.0e-9*sqrt(phys.T/300);
            end
        case 'F'
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 2)
                as_drv_cv.noise.Vn_opa = 26.4e-9*sqrt(phys.T/300);
            end
            if (asic_config.otp.ana.gyr_trm_cvd_gm == 3)
                as_drv_cv.noise.Vn_opa = 35.0e-9*sqrt(phys.T/300); 
            end
    end

end
