function [qc, gyr_trm_cap_prog_qc_ch1, gyr_trm_cap_prog_qc_ch2, gyr_trm_cap_prog_qc_ch3] = init_qc(varargin)
    asic_struct = varargin{1};
    asic_config = varargin{2};
    if length(varargin) > 2
        gyro_struct = varargin{3};
    end
    qc = struct;
    if length(varargin) > 2
        gyro_sens_cv_out_slope_ch1 = ...
            gyro_struct.stat.mc.Sdd_C_Fdps * (asic_config.inf.ref.cm.vout - asic_config.inf.sense.cv.surface_charge) / asic_config.inf.sense.cv.C.gyr_cap_cv_s ...
            * asic_config.inf.drive.cv.drive_ampl / gyro_struct.stat.mc.xd ;
        gyro_sens_cv_out_slope_ch2 = ...
            gyro_struct.stat.mc.Sdd_S_Fdps * (asic_config.inf.ref.cm.vout - asic_config.inf.sense.cv.surface_charge) / asic_config.inf.sense.cv.S.gyr_cap_cv_s ...
            * asic_config.inf.drive.cv.drive_ampl / gyro_struct.stat.mc.xd ;
        gyro_sens_cv_out_slope_ch3 = ...
            gyro_struct.stat.mc.Sdd_Z_Fdps * (asic_config.inf.ref.cm.vout - asic_config.inf.sense.cv.surface_charge) / asic_config.inf.sense.cv.Z.gyr_cap_cv_s ...
            * asic_config.inf.drive.cv.drive_ampl / gyro_struct.stat.mc.xd ;

        gyr_trm_cap_prog_qc_ch1 = round(sign(gyro_struct.stat.mc.fsplit_C)*gyro_struct.stat.mc.quad_C * gyro_struct.stat.mc.Sdd_C_Fdps * asic_config.inf.ref.cm.vout / asic_config.inf.drive.agc.out / asic_struct.sense.fe.cv.cap_qc(2) *(1 + asic_config.otp.ana.gyr_trm_drv_4um_en));
        gyr_trm_cap_prog_qc_ch2 = round(sign(gyro_struct.stat.mc.fsplit_S)*gyro_struct.stat.mc.quad_S * gyro_struct.stat.mc.Sdd_S_Fdps * asic_config.inf.ref.cm.vout / asic_config.inf.drive.agc.out / asic_struct.sense.fe.cv.cap_qc(2) *(1 + asic_config.otp.ana.gyr_trm_drv_4um_en));
        gyr_trm_cap_prog_qc_ch3 = round(sign(gyro_struct.stat.mc.fsplit_Z)*gyro_struct.stat.mc.quad_Z * gyro_struct.stat.mc.Sdd_Z_Fdps * asic_config.inf.ref.cm.vout / asic_config.inf.drive.agc.out / asic_struct.sense.fe.cv.cap_qc(2) *(1 + asic_config.otp.ana.gyr_trm_drv_4um_en));
        
        qc.C.gyr_cap_qc = sum(asic_struct.sense.fe.cv.cap_qc.*dec2binvec(abs(gyr_trm_cap_prog_qc_ch1), 9, 0));
        qc.S.gyr_cap_qc = sum(asic_struct.sense.fe.cv.cap_qc.*dec2binvec(abs(gyr_trm_cap_prog_qc_ch2), 9, 0));
        qc.Z.gyr_cap_qc = sum(asic_struct.sense.fe.cv.cap_qc.*dec2binvec(abs(gyr_trm_cap_prog_qc_ch3), 9, 0));
        
        % calculate back which quad is really compensated
        qc.C.quad_comp = gyr_trm_cap_prog_qc_ch1 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.C.gyr_cap_cv_s;
        qc.S.quad_comp = gyr_trm_cap_prog_qc_ch2 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.S.gyr_cap_cv_s;
        qc.Z.quad_comp = gyr_trm_cap_prog_qc_ch3 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.Z.gyr_cap_cv_s;
        
        qc.C.quad_comp_dps = gyr_trm_cap_prog_qc_ch1 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.C.gyr_cap_cv_s / gyro_sens_cv_out_slope_ch1;
        qc.S.quad_comp_dps = gyr_trm_cap_prog_qc_ch2 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.S.gyr_cap_cv_s / gyro_sens_cv_out_slope_ch2;
        qc.Z.quad_comp_dps = gyr_trm_cap_prog_qc_ch3 * asic_struct.sense.fe.cv.cap_qc(2) * asic_config.inf.drive.agc.out / asic_config.inf.sense.cv.Z.gyr_cap_cv_s / gyro_sens_cv_out_slope_ch3;
       
        % Residual of quadrature after quad compensation [dps]
        qc.C.quad_res = abs(gyro_struct.stat.mc.quad_C) - abs(qc.C.quad_comp_dps) + gyro_struct.stat.mc.quad_drift_C;
        qc.S.quad_res = abs(gyro_struct.stat.mc.quad_S) - abs(qc.S.quad_comp_dps) + gyro_struct.stat.mc.quad_drift_S;
        qc.Z.quad_res = abs(gyro_struct.stat.mc.quad_Z) - abs(qc.Z.quad_comp_dps) + gyro_struct.stat.mc.quad_drift_Z;

    end
end
