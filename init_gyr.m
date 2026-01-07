function init_gyr(obj)
    %fixdp_lut
    obj.sense.gyr_runit_rppoly_scale = 1;
    obj.cm.avdd = 1.5;       % [V] analog supply voltage (1.5V)
    obj.cm.vref = 1.210;     % [V] analog reference voltage
    
    obj.cm.vgm  = ...
        [  0, 8.9286;  1, 9.0435;  2, 9.1618;  3, 9.2836;  4, 9.4091;  5, 9.5385;  6, 9.6719; ...
           7, 9.8095;  8, 9.9516;  9, 10.0984; 10, 10.2500; 11, 10.4068; 12, 10.5690; 13, 10.7368];
    
    
    %% drive vco parameters
    obj.drive.vco.i_curr = 5e-6;          % [A] vco current
    obj.drive.vco.cap_par = 30.0e-015;    % full parasitic capacity of vco stage provided by neg7rt
    %% 
    
    %% TODO: amplifier gains and corrections 10^(g0/20)*(1+3*gb1*u^6+gb1*u^8)
    obj.drive.fe.cv.g0  = 70;        % TODO: double-check
    obj.drive.fe.cv.gb1 = -0.007388; % TODO: double-check (normally extracted from analog designers)
    
    obj.drive.fe.agc.g0  = 70;
    obj.drive.fe.agc.gb1 = -0.007388;
    
    obj.drive.fe.pll.g0  = 70;
    obj.drive.fe.pll.gb1 = -0.007388;
    
    obj.sense.fe.cv.g0  = 108.8;
    obj.sense.fe.cv.gb1 = -0.004594;
    
    
    obj.sense.fe.rate.g0  = 108.8;
    obj.sense.fe.rate.gb1 = -0.004594;
    
    obj.sense.fe.quad.g0  = 108.8;
    obj.sense.fe.quad.gb1 = -0.004594;
    
    
    obj.sense.fe.temp.g0  = 108.8;
    obj.sense.fe.temp.gb1 = -0.004594;
    
    %% TODO: sense channel, CV-Converter offset and DCS feedback ratio
    
    obj.sense.fe.cv.Voffs = 0.2;
    obj.sense.fe.cv.alfa_fb = 0.9;
    
    %% parasitic capacitor between drive and sense channels
    % TODO: Extract the correct values
    
    obj.drive.fe.cv.cp_csp_aa = 0.00e-14;   % [F]
    obj.drive.fe.cv.cp_csp_ai = 0.00e-15;   % [F]
    obj.drive.fe.cv.cp_csn_aa = 0.00e-14;   % [F]
    obj.drive.fe.cv.cp_csn_ai = 0.00e-14;   % [F]
    
    obj.sense.fe.x.cp_csp_aa = 0.00e-15;    % [F]
    obj.sense.fe.x.cp_csp_ai = 0.00e-15;    % [F]
    obj.sense.fe.x.cp_csn_aa = 0.00e-15;    % [F]
    obj.sense.fe.x.cp_csn_ai = 0.00e-15;    % [F]
    
    obj.sense.fe.y.cp_csp_aa = 0.00e-15;    % [F]
    obj.sense.fe.y.cp_csp_ai = 0.00e-15;    % [F]
    obj.sense.fe.y.cp_csn_aa = 0.00e-15;    % [F]
    obj.sense.fe.y.cp_csn_ai = 0.00e-15;    % [F]
    
    obj.sense.fe.z.cp_csp_aa = 0.00e-15;    % [F]
    obj.sense.fe.z.cp_csp_ai = 0.00e-15;    % [F]
    obj.sense.fe.z.cp_csn_aa = 0.00e-15;    % [F]
    obj.sense.fe.z.cp_csn_ai = 0.00e-15;    % [F]
    
    obj.sense.fe.cv.cp_csp_aa = [obj.sense.fe.x.cp_csp_aa;obj.sense.fe.y.cp_csp_aa;obj.sense.fe.z.cp_csp_aa];
    obj.sense.fe.cv.cp_csp_ai = [obj.sense.fe.x.cp_csp_ai;obj.sense.fe.y.cp_csp_ai;obj.sense.fe.z.cp_csp_ai];
    obj.sense.fe.cv.cp_csn_aa = [obj.sense.fe.x.cp_csn_aa;obj.sense.fe.y.cp_csn_aa;obj.sense.fe.z.cp_csn_aa];
    obj.sense.fe.cv.cp_csn_ai = [obj.sense.fe.x.cp_csn_ai;obj.sense.fe.y.cp_csn_ai;obj.sense.fe.z.cp_csn_ai];
    
    %% resistor and capacitor values in the drive loop (don't change!!!)
    
    % capacitor values for c/v gain
    obj.drive.fe.cv.cap_cv         = struct;
    obj.drive.fe.cv.cap_cv.c_fixed = 63.275e-15*16; % when DRV_4UM_EN=0
    obj.drive.fe.cv.cap_cv.c_bank0 = 63.275e-15*8;  % when DRV_4UM_EN=1
    obj.drive.fe.cv.cap_cv.c_bank1 = 63.275e-15*1;  
    obj.drive.fe.cv.cap_cv.c_bank2 = 63.275e-15*2;
    obj.drive.fe.cv.cap_cv.c_bank3 = 63.275e-15*4;
    obj.drive.fe.cv.cap_cv.c_bank4 = 63.275e-15*8;
    obj.drive.fe.cv.cap_cv.c_bank5 = 63.275e-15*16; % when DRV_4UM_EN=0
    obj.drive.fe.cv.cap_cv.c_bank6 = 0.0e-15;       % when DRV_4UM_EN=1
    obj.drive.fe.cv.cap_cv         = cell2mat(struct2cell( obj.drive.fe.cv.cap_cv ));
    
     % resistor and capacitor values for gain in agc integrator
    obj.drive.fe.agc.cap_int         = struct;
    obj.drive.fe.agc.cap_int.c_fixed = 480.715e-15*10;
    obj.drive.fe.agc.cap_int.c_bank0 = 480.715e-15*10;
    obj.drive.fe.agc.cap_int.c_bank1 = 480.715e-15*10;
    obj.drive.fe.agc.cap_int.c_bank2 = 0.0e-12;
    obj.drive.fe.agc.cap_int.c_bank3 = 0.0e-12;
    obj.drive.fe.agc.cap_int.c_bank4 = 0.0e-12;
    obj.drive.fe.agc.cap_int         = cell2mat(struct2cell( obj.drive.fe.agc.cap_int ));
    
    obj.drive.fe.agc.res_sig         = struct;
    obj.drive.fe.agc.res_sig.r_fixed = 1170.02e3+390.006e3;
    obj.drive.fe.agc.res_sig.r_bank0 = 24.375e3*1;
    obj.drive.fe.agc.res_sig.r_bank1 = 24.375e3*2;
    obj.drive.fe.agc.res_sig.r_bank2 = 24.375e3*4;
    obj.drive.fe.agc.res_sig.r_bank3 = 24.375e3*8;
    obj.drive.fe.agc.res_sig.r_bank4 = 24.375e3*16;
    obj.drive.fe.agc.res_sig.r_bank5 = 24.375e3*32;
    obj.drive.fe.agc.res_sig         = cell2mat(struct2cell(obj.drive.fe.agc.res_sig));
    
    obj.drive.fe.agc.res_ref         = struct;
    obj.drive.fe.agc.res_ref.r_fixed = 1247.5e3*2;
    obj.drive.fe.agc.res_ref.r_bank0 = 38.984e3;
    obj.drive.fe.agc.res_ref.r_bank1 = 38.984e3*2;
    obj.drive.fe.agc.res_ref.r_bank2 = 38.984e3*4;
    obj.drive.fe.agc.res_ref.r_bank3 = 38.984e3*8;
    obj.drive.fe.agc.res_ref.r_bank4 = 38.984e3*16;
    obj.drive.fe.agc.res_ref.r_bank5 = 38.984e3*32;
    obj.drive.fe.agc.res_ref         = cell2mat(struct2cell( obj.drive.fe.agc.res_ref ));
    
    % resistor and capacitor values for gain in phase integrator 
    obj.drive.fe.pll.cap_int         = struct;
    obj.drive.fe.pll.cap_int.c_fixed = 480.715e-15*14; % +480.715e-15*15 when LOW_BW_EN=1
    obj.drive.fe.pll.cap_int.c_bank0 = 480.715e-15/2;
    obj.drive.fe.pll.cap_int.c_bank1 = 480.715e-15;
    obj.drive.fe.pll.cap_int.c_bank2 = 480.715e-15*2;
    obj.drive.fe.pll.cap_int.c_bank3 = 480.715e-15*4;
    obj.drive.fe.pll.cap_int.c_bank4 = 480.715e-15*8;
    obj.drive.fe.pll.cap_int         = cell2mat(struct2cell( obj.drive.fe.pll.cap_int ));
    
    obj.drive.fe.pll.res_sig         = struct;
    obj.drive.fe.pll.res_sig.r_fixed = 1912.830e3;
    obj.drive.fe.pll.res_sig.r_bank0 = 30.3565e3;
    obj.drive.fe.pll.res_sig.r_bank1 = 30.3565e3*2;
    obj.drive.fe.pll.res_sig.r_bank2 = 30.3565e3*4;
    obj.drive.fe.pll.res_sig.r_bank3 = 30.3565e3*8;
    obj.drive.fe.pll.res_sig.r_bank4 = 30.3565e3*16;
    obj.drive.fe.pll.res_sig.r_bank5 = 30.3565e3*32;
    obj.drive.fe.pll.res_sig         = cell2mat(struct2cell( obj.drive.fe.pll.res_sig ));
    
    %not used, but needed
    obj.drive.fe.pll.res_ref         = struct;
    obj.drive.fe.pll.res_ref.r_fixed = 1366.3e3;
    obj.drive.fe.pll.res_ref.r_bank0 = 21.683e3;
    obj.drive.fe.pll.res_ref.r_bank1 = 21.683e3*2;
    obj.drive.fe.pll.res_ref.r_bank2 = 21.683e3*4;
    obj.drive.fe.pll.res_ref.r_bank3 = 21.683e3*8;
    obj.drive.fe.pll.res_ref.r_bank4 = 21.683e3*16;
    obj.drive.fe.pll.res_ref.r_bank5 = 21.683e3*32;
    obj.drive.fe.pll.res_ref         = cell2mat(struct2cell( obj.drive.fe.pll.res_ref ));
    
    % resistor and capacitor values for gain in vco
    obj.drive.vco.cap_vco         = struct;
    obj.drive.vco.cap_vco.c_fixed = (37.1083e-15*5+30.0e-15)*2;
    obj.drive.vco.cap_vco.c_bank0 = 37.1083e-15/2;
    obj.drive.vco.cap_vco.c_bank1 = 37.1083e-15;
    obj.drive.vco.cap_vco.c_bank2 = 37.1083e-15*2;
    obj.drive.vco.cap_vco.c_bank3 = 37.1083e-15*4;
    obj.drive.vco.cap_vco.c_bank4 = 37.1083e-15*8;
    obj.drive.vco.cap_vco         = cell2mat(struct2cell( obj.drive.vco.cap_vco ));
    
    obj.drive.vco.res_sig         = struct;
    obj.drive.vco.res_sig.r_fixed = 1011.37e3+252.843e3;
    obj.drive.vco.res_sig.r_bank0 = 31.6053e3;
    obj.drive.vco.res_sig.r_bank1 = 31.6053e3*2;
    obj.drive.vco.res_sig.r_bank2 = 31.6053e3*4;
    obj.drive.vco.res_sig.r_bank3 = 31.6053e3*8;
    obj.drive.vco.res_sig.r_bank4 = 31.6053e3*16;
    obj.drive.vco.res_sig.r_bank5 = 1011.37e3; % ~31.6053e3*32 
    obj.drive.vco.res_sig         = cell2mat(struct2cell( obj.drive.vco.res_sig ));
    
    obj.drive.vco.res_ref         = struct;
    obj.drive.vco.res_ref.r_fixed = 210.097e3/3;
    obj.drive.vco.res_ref.r_bank0 = 13.131e3/12;
    obj.drive.vco.res_ref.r_bank1 = 13.131e3/6;
    obj.drive.vco.res_ref.r_bank2 = 13.131e3*2/6;
    obj.drive.vco.res_ref.r_bank3 = 13.131e3*4/6;
    obj.drive.vco.res_ref.r_bank4 = 13.131e3*8/6;
    obj.drive.vco.res_ref.r_bank5 = 13.131e3*16/6;
    obj.drive.vco.res_ref         = cell2mat(struct2cell( obj.drive.vco.res_ref ));
    
    
    %% resistor and capacitor values in the sense part
    
    % capacitor values for c/v gain;
    obj.sense.fe.cv.cap_cv         = struct;
    obj.sense.fe.cv.cap_cv.c_fixed = 47.7309e-15;    % when FS_2X_EN=0
    obj.sense.fe.cv.cap_cv.c_bank0 = 47.7309e-15*2;  % when FS_2X_EN=1
    obj.sense.fe.cv.cap_cv.c_bank1 =  8.0336e-15;    % when FS_2X_EN=0
    obj.sense.fe.cv.cap_cv.c_bank2 =  8.0336e-15*2; 
    obj.sense.fe.cv.cap_cv.c_bank3 =  8.0336e-15*4;
    obj.sense.fe.cv.cap_cv.c_bank4 =  8.0336e-15*8;
    obj.sense.fe.cv.cap_cv.c_bank5 =  8.0336e-15*16;
    obj.sense.fe.cv.cap_cv.c_bank6 =  8.0336e-15*32; % when FS_2X_EN=1
    obj.sense.fe.cv.cap_cv         = cell2mat(struct2cell( obj.sense.fe.cv.cap_cv ));
    
    %% all zero is disabled, 104 dec for 24kdps
    obj.sense.fe.cv.cap_qc         = struct;
    obj.sense.fe.cv.cap_qc.c_fixed = 0.0;
    obj.sense.fe.cv.cap_qc.c_bank0 = 4.0168e-15;
    obj.sense.fe.cv.cap_qc.c_bank1 = 4.0168e-15*2;
    obj.sense.fe.cv.cap_qc.c_bank2 = 4.0168e-15*4; 
    obj.sense.fe.cv.cap_qc.c_bank3 = 4.0168e-15*8;
    obj.sense.fe.cv.cap_qc.c_bank4 = 4.0168e-15*16;
    obj.sense.fe.cv.cap_qc.c_bank5 = 4.0168e-15*32;
    obj.sense.fe.cv.cap_qc.c_bank6 = 4.0168e-15*64;
    obj.sense.fe.cv.cap_qc.c_bank7 = 4.0168e-15*128;
    obj.sense.fe.cv.cap_qc.c_bank8 = 0.0; %sign
    obj.sense.fe.cv.cap_qc         = cell2mat(struct2cell( obj.sense.fe.cv.cap_qc ));
    
    
    %% resistor and capacitor values for gain in rate-path;
    
    obj.sense.fe.rate.res_sig1         = struct;
    obj.sense.fe.rate.res_sig1.r_fixed = 1912.837e3;
    obj.sense.fe.rate.res_sig1.r_bank0 = 30.3565e3;
    obj.sense.fe.rate.res_sig1.r_bank1 = 30.3565e3*2;
    obj.sense.fe.rate.res_sig1.r_bank2 = 30.3565e3*4;
    obj.sense.fe.rate.res_sig1.r_bank3 = 30.3565e3*8;
    obj.sense.fe.rate.res_sig1.r_bank4 = 30.3565e3*16;
    obj.sense.fe.rate.res_sig1.r_bank5 = 30.3565e3*32;
    obj.sense.fe.rate.res_sig1         = cell2mat(struct2cell( obj.sense.fe.rate.res_sig1 ));

    obj.sense.fe.rate.res_sig2         = struct;
    obj.sense.fe.rate.res_sig2.r_fixed = 1366.31e3;
    obj.sense.fe.rate.res_sig2.r_bank0 = 21.6832e3;
    obj.sense.fe.rate.res_sig2.r_bank1 = 21.6832e3*2;
    obj.sense.fe.rate.res_sig2.r_bank2 = 21.6832e3*4;
    obj.sense.fe.rate.res_sig2.r_bank3 = 21.6832e3*8;
    obj.sense.fe.rate.res_sig2.r_bank4 = 21.6832e3*16;
    obj.sense.fe.rate.res_sig2.r_bank5 = 21.6832e3*32;
    obj.sense.fe.rate.res_sig2         = cell2mat(struct2cell( obj.sense.fe.rate.res_sig2 ));

    % res_offset = 203.209e3
    
    obj.sense.fe.rate.res_fb         = struct;
    obj.sense.fe.rate.res_fb.r_fixed = 1366.31e3;
    obj.sense.fe.rate.res_fb.r_bank0 = 21.6832e3;
    obj.sense.fe.rate.res_fb.r_bank1 = 21.6832e3*2;
    obj.sense.fe.rate.res_fb.r_bank2 = 21.6832e3*4;
    obj.sense.fe.rate.res_fb.r_bank3 = 21.6832e3*8;
    obj.sense.fe.rate.res_fb.r_bank4 = 21.6832e3*16;
    obj.sense.fe.rate.res_fb.r_bank5 = 21.6832e3*32;
    obj.sense.fe.rate.res_fb         = cell2mat(struct2cell( obj.sense.fe.rate.res_fb ));
    
    % TODO: update res_quadint
    obj.sense.fe.rate.res_quadint         = struct;
    obj.sense.fe.rate.res_quadint.r_fixed = 812.837e3;
    obj.sense.fe.rate.res_quadint.r_bank0 = 12.897e3;
    obj.sense.fe.rate.res_quadint.r_bank1 = 25.809e3;
    obj.sense.fe.rate.res_quadint.r_bank2 = 51.608e3;
    obj.sense.fe.rate.res_quadint.r_bank3 = 103.207e3;
    obj.sense.fe.rate.res_quadint.r_bank4 = 206.414e3;
    obj.sense.fe.rate.res_quadint.r_bank5 = 412.828e3;
    obj.sense.fe.rate.res_quadint         = cell2mat(struct2cell( obj.sense.fe.rate.res_quadint ));
    
    % due to  technology
    obj.sense.fe.rate.cap_int         = struct;
    obj.sense.fe.rate.cap_int.c_fixed = 350.0e-15*29;
    obj.sense.fe.rate.cap_int.c_bank0 = 350.0e-15;
    obj.sense.fe.rate.cap_int.c_bank1 = 350.0e-15*2;
    obj.sense.fe.rate.cap_int.c_bank2 = 350.0e-15*4;
    obj.sense.fe.rate.cap_int.c_bank3 = 350.0e-15*8;
    obj.sense.fe.rate.cap_int.c_bank4 = 0.0;
    obj.sense.fe.rate.cap_int         = cell2mat(struct2cell( obj.sense.fe.rate.cap_int ));
    
    obj.sense.be.rate.quant_gain = 4.5;
    obj.sense.be.rate.quant_vref = 1.210;
    
    %% idac reference trim
    obj.cm.idac.res_ref          = struct;
    obj.cm.idac.res_ref.r_fixed  = 45.5437e3*2 + 182.175e3;
    obj.cm.idac.res_ref.r_bank0  = 4.33664e3;
    obj.cm.idac.res_ref.r_bank1  = 8.67387e3;
    obj.cm.idac.res_ref.r_bank2  = 17.3483e3;
    obj.cm.idac.res_ref.r_bank3  = 34.6973e3;
    obj.cm.idac.res_ref.r_bank4  = 69.3946e3;
    obj.cm.idac.res_ref.r_bank5  = 138.789e3;
    obj.cm.idac.res_ref          = cell2mat(struct2cell( obj.cm.idac.res_ref ));
  
    % TODO: Move to asic_trim ff and fb
    obj.sense.be.rate.a1  = 6/6;
    obj.sense.be.rate.ff1 = 1/2;
    obj.sense.be.rate.fb1 = 2/6;
    
    obj.sense.be.rate.a2  = 4/6;
    obj.sense.be.rate.fb2 = 2/8;
    obj.sense.be.rate.fb3 = 1/6;
    
    obj.sense.be.rate.d1  = 6/6; % TODO: totally move from here as well as from be_rate model
    obj.sense.be.rate.d2  = 5/2;
    
    obj.sense.be.rate.s1 = 0;
    
    obj.sense.be.rate.c1 = 1/8;
    obj.sense.be.rate.c2 = 1/2;
    
    %TODO: not important for now
    obj.sense.be.rate.dith_en    = 1;
    obj.sense.be.rate.dith_power = 0.0092;
    
    %% resistor and capacitor values for gain in quad-path; 
    obj.sense.fe.quad.res_sig         = struct;
    obj.sense.fe.quad.res_sig.r_fixed = 1912.837e3;
    obj.sense.fe.quad.res_sig.r_bank0 = 30.3565e3;
    obj.sense.fe.quad.res_sig.r_bank1 = 30.3565e3*2;
    obj.sense.fe.quad.res_sig.r_bank2 = 30.3565e3*4;
    obj.sense.fe.quad.res_sig.r_bank3 = 30.3565e3*8;
    obj.sense.fe.quad.res_sig.r_bank4 = 30.3565e3*16;
    obj.sense.fe.quad.res_sig.r_bank5 = 30.3565e3*32;
    obj.sense.fe.quad.res_sig         = cell2mat(struct2cell( obj.sense.fe.quad.res_sig ));
    
    obj.sense.fe.quad.res_fb         = struct;
    obj.sense.fe.quad.res_fb.r_fixed = 1366.31e3;
    obj.sense.fe.quad.res_fb.r_bank0 = 21.6832e3;
    obj.sense.fe.quad.res_fb.r_bank1 = 21.6832e3*2;
    obj.sense.fe.quad.res_fb.r_bank2 = 21.6832e3*4;
    obj.sense.fe.quad.res_fb.r_bank3 = 21.6832e3*8;
    obj.sense.fe.quad.res_fb.r_bank4 = 21.6832e3*16;
    obj.sense.fe.quad.res_fb.r_bank5 = 21.6832e3*32;
    obj.sense.fe.quad.res_fb         = cell2mat(struct2cell( obj.sense.fe.quad.res_fb ));
    
    % TODO: update res_quadint
    obj.sense.fe.quad.res_quadint         = struct;
    obj.sense.fe.quad.res_quadint.r_fixed = 812.837e3;
    obj.sense.fe.quad.res_quadint.r_bank0 = 12.897e3;
    obj.sense.fe.quad.res_quadint.r_bank1 = 25.809e3;
    obj.sense.fe.quad.res_quadint.r_bank2 = 51.608e3;
    obj.sense.fe.quad.res_quadint.r_bank3 = 103.207e3;
    obj.sense.fe.quad.res_quadint.r_bank4 = 206.414e3;
    obj.sense.fe.quad.res_quadint.r_bank5 = 412.828e3;
    obj.sense.fe.quad.res_quadint         = cell2mat(struct2cell( obj.sense.fe.quad.res_quadint ));
    
    obj.sense.fe.quad.cap_int = struct;
    obj.sense.fe.quad.cap_int.c_fixed = 350.0e-15*29;
    obj.sense.fe.quad.cap_int.c_bank0 = 350.0e-15*2;
    obj.sense.fe.quad.cap_int.c_bank1 = 350.0e-15*4;
    obj.sense.fe.quad.cap_int.c_bank2 = 350.0e-15*8;
    obj.sense.fe.quad.cap_int.c_bank3 = 350.0e-15*16;
    obj.sense.fe.quad.cap_int.c_bank4 = 0.0;
    obj.sense.fe.quad.cap_int         = cell2mat(struct2cell( obj.sense.fe.quad.cap_int ));
    
    obj.sense.be.quad.quant_gain = 4.5;
    obj.sense.be.quad.quant_vref = 1.210;
    
    obj.sense.be.quad.a1  = 6/6;
    obj.sense.be.quad.ff1 = 1/2;
    obj.sense.be.quad.fb1 = 2/6;
    
    obj.sense.be.quad.a2  = 4/6;
    obj.sense.be.quad.fb2 = 2/8;
    obj.sense.be.quad.fb3 = 1/6;
    
    obj.sense.be.quad.d1  = 6/6;
    obj.sense.be.quad.d2  = 5/2;
    
    obj.sense.be.quad.s1 = 0;
    
    obj.sense.be.quad.c1 = 1/8;
    obj.sense.be.quad.c2 = 1/2;
    
    obj.sense.be.quad.dith_en    = 0;
    obj.sense.be.quad.dith_power = 0.0092;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
