function otp = init_gyr_otp()
    %% ana reg default config bits
    otp.ana.gyr_trm_cap_prog_cv_d         = hex2dec('09');
    otp.ana.gyr_trm_adc_int_dem_sel       = hex2dec('0');
    otp.ana.gyr_trm_adc_single_bit_q_sel  = hex2dec('0');
    otp.ana.gyr_trm_adc_single_bit_r_sel  = hex2dec('0');
    otp.ana.gyr_trm_cvd_gm                = hex2dec('2'); %fi(2,0,2,0); %
    otp.ana.gyr_trm_cvs_gm_ch1            = hex2dec('2'); %fi(2,0,2,0); %
    otp.ana.gyr_trm_cvs_gm_ch2            = hex2dec('2'); %fi(2,0,2,0); %
    otp.ana.gyr_trm_cvs_gm_ch3            = hex2dec('2'); %fi(2,0,2,0); %
    otp.ana.gyr_trm_res_prog_bg           = hex2dec('22');%fi(34,0,8,0);%
    otp.ana.gyr_trm_rint_ib               = hex2dec('0');
    otp.ana.gyr_trm_pll_lowbw_en          = hex2dec('0');
    otp.ana.gyr_trm_cap_prog_a            = hex2dec('0');
    otp.ana.gyr_trm_drv_4um_en            = hex2dec('0');
    otp.ana.gyr_trm_lp_dpint_en           = hex2dec('1');
    otp.ana.gyr_trm_phbe_lp_en            = hex2dec('1');
    otp.ana.gyr_trm_invert_drive          = hex2dec('0');
    otp.ana.gyr_trm_cap_prog_r            = hex2dec('10');
    otp.ana.gyr_trm_cap_prog_p            = hex2dec('10');
    otp.ana.gyr_trm_cap_prog_q            = hex2dec('10');
    otp.ana.gyr_trm_dith_gain_sd_q        = hex2dec('0');
    otp.ana.gyr_trm_dith_gain_sd_r        = hex2dec('0');
    otp.ana.gyr_trm_daint_chop_dis        = hex2dec('0');
    otp.ana.gyr_trm_clk_sel_bp            = hex2dec('0');
    otp.ana.gyr_trm_sns_be_lp_en          = hex2dec('1');
    otp.ana.gyr_trm_cvs_fs_2x_en          = hex2dec('0');
    otp.ana.gyr_test_cvs_highfs_en        = hex2dec('0');
    otp.ana.gyr_trm_spare                 = hex2dec('0');
    otp.ana.gyr_trm_res_prog_ldovcm       = hex2dec('8');
    otp.ana.gyr_trm_res_prog_ldo_a        = hex2dec('13');

    %% fcu reg default config bits
    otp.fcu.c_odr_trim                    = hex2dec('0');
    otp.fcu.high_prec_mode                = hex2dec('0');
    otp.fcu.gyr_config_mash_2_1_en        = hex2dec('0');
    otp.fcu.gyr_odr                       = hex2dec('D');
    otp.fcu.gyr_range                     = hex2dec('5');
    otp.fcu.gyr_bw                        = hex2dec('2');
    otp.fcu.gyr_avg_num                   = hex2dec('2');
    otp.fcu.gyr_mode                      = hex2dec('7');
    otp.fcu.gyr_filter_type               = hex2dec('1');
    otp.fcu.gyr_filter_latency            = hex2dec('0');
    otp.fcu.gyr_notch_byp                 = hex2dec('0');
    otp.fcu.gyr_test_high_odr_en          = hex2dec('0');
    otp.fcu.gyr_st_force_dis              = hex2dec('0');
    otp.fcu.gyr_test_quad_mode            = hex2dec('0');

    %% fcu reg drive duty-cycle bits
    otp.fcu.gyr_drv_margin                = hex2dec('20');
    
    %% gyr_reg default config bits
    otp.gyr_afe.gyr_notch_grp_comp_en     = hex2dec('1');
    otp.gyr_afe.gyr_afe_wait              = hex2dec('4');
    otp.gyr_afe.gyr_grp_comp_length       = hex2dec('0');
    otp.gyr_afe.gyr_average_mode          = hex2dec('0');
    otp.gyr_afe.gyr_mash_bypass           = hex2dec('0');
    otp.gyr_afe.gyr_mash_dly_bypass       = hex2dec('0');
    otp.gyr_afe.gyr_mash_dgain_d1         = hex2dec('0');
    otp.gyr_afe.gyr_mash_dgain_d2         = hex2dec('5');
    otp.gyr_afe.gyr_test_tosk_en          = hex2dec('0');
    otp.gyr_afe.gyr_test_tosk_qp          = hex2dec('0');
    otp.gyr_afe.gyr_test_tosk_qn          = hex2dec('0');
    otp.gyr_afe.gyr_test_iqmod_en         = hex2dec('0');
    otp.gyr_afe.gyr_test_iqmod_sin        = hex2dec('0');
    otp.gyr_afe.gyr_test_iqmod_cos        = hex2dec('0');
    otp.gyr_afe.gyr_test_ac2qpqn_en       = hex2dec('0');
    otp.gyr_afe.gyr_test_ac2qpqn_sel_rate = hex2dec('0');
    otp.gyr_afe.gyr_test_quad_sel         = hex2dec('0');
    otp.gyr_afe.gyr_test_quad2rate_sel    = hex2dec('0');
    otp.gyr_afe.gyr_test_quad_avg_bypass  = hex2dec('0');
    otp.gyr_afe.gyr_test_quad_out_mute    = hex2dec('0');
    otp.gyr_afe.gyr_test_rdem_sel         = hex2dec('0');
    otp.gyr_afe.gyr_test_qdem_sel         = hex2dec('0');

    %% drv_st_reg default config bits
    otp.drv_st.gyr_mar_aout_diff_en       = hex2dec('0');
    otp.drv_st.gyr_mar_aout_diff          = hex2dec('0');
    otp.drv_st.gyr_mar_pout_diff          = hex2dec('0');
    otp.drv_st.gyr_mar_drv_mag            = hex2dec('0');
    otp.drv_st.gyr_filt_bite_res          = hex2dec('0');
    otp.drv_st.gyr_bite_startup_dis       = hex2dec('0');
  
    %% VCO
    otp.ana.gyr_trm_res_prog_vco          = hex2dec('20');    
    otp.ana.gyr_trm_res_prog_i            = hex2dec('10');
    otp.ana.gyr_trm_res_prog_i2           = hex2dec('0');
end
