clc
set(0,...
    'defaultFigureColor','w',...
    'DefaultFigureWindowStyle','docked', ...
    'defaultTextInterpreter','tex', ...
    'defaultAxesFontSize',14, ...
    'defaultLegendInterpreter','tex', ...
    'defaultLineLineWidth', 1);

[config, logging] = init_sim();
config.sys.temp.Value = 300;
config.sys.rate_x.Value = 10.0;
config.sys.rate_y.Value = 100.0;
config.sys.rate_z.Value = 1000.0;

config.sys.gyr.mode = 'LNM';
config.sys.filt_type = 'IIR';
config.sys.has_notch = [0; 0];
config.sys.odr_sim = 'SINGLE';
config.aliasing_en = 0;

switch config.sys.filt_type
    case {'BYPASS'}
        config.sys.odr_list = [6400; 3200; 1600; 800];
        config.sys.idr_list = [6400; 3200; 1600; 800];
        config.sys.bw_list =  [];
        config.sys.rolloff_list = [];
    case {'IIR', 'IIR_W_ZERO'}
        switch config.sys.odr_sim
            case {'FULL'}
                config.sys.odr_list = [
                    %6400; 6400; 6400; 6400; ...
                    3200; ...
                    ];
                config.sys.idr_list = [
                    %6400; 6400; 6400; 6400; ...
                    3200; ...
                    ];
                config.sys.bw_list = [
                    %1600;  800;  400;  200; ...
                    100; ...
                    ];
                config.sys.rolloff_list = [20; 40; 60];
            case {'SINGLE'}
                config.sys.odr_list = 6400;
                config.sys.idr_list = 6400;
                config.sys.bw_list = 800;
                config.sys.rolloff_list = 60;
        end
        config.spec.param.latency_ms = 30.7/360/20*1000;
        config.spec.param.latency_freq = 20;
        config.spec.param.mag_3dB = -3;
        
    case {'AVG'}
        config.sys.odr_list = 6400;
        config.sys.bw_list =  [];
        config.sys.rolloff_list = [];
end
mc10000 = [
    "warp_g210_p254_10000runs_240614_142034_IO.txt", ...
    "warp_g211_p122_10000runs_240619_113259_IO.txt", ...
    "warp_g210_p312_10000runs_240614_160004_IO.txt", ...
    "warp_g210_p507_10000runs_240614_123251_IO.txt", ...
    "warp_g210_p249_10000runs_240614_102957_IO.txt", ...
    "warp_g210_p408_10000runs_240617_113900_IO.txt", ...
    "warp_g214_p288_wMSA_10000runs_241024_213300_IO.csv"
    ];

mems_par = 'PEX_Gyro_C10.xlsx';
mems_sheet = 'Gyro_C10_20241030_1';
system0 = imu(1,220,config, mc10000(7), mems_par, mems_sheet, 'HPM');
simulink_data = system0.gyro.gyro.simulink_data;
clear mc10000;
% isnoisy = 0;
% bai430aa_noisy = Simulink.VariantExpression('isnoisy == 1');
% bai430aa_noiseless = Simulink.VariantExpression('isnoisy == 0');

% mex -v bmi430_prj/vco_fcn.c
% mex -v bmi430_prj/analysis/noise/delsig/simulateDSM_mod.c
% mex -v bmi430_prj/analysis/noise/delsig/simulateDSM.c
mex -v bmi430_prj/vco_fcn.c
mex -v bmi430_prj/analysis/noise/delsig/simulateDSM_mod.c
mex -v bmi430_prj/analysis/noise/delsig/simulateDSM.c
