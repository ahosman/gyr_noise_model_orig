function [Sdd_C, Sdd_S, Sdd_Z, stot_C, stot_S, stot_Z, gyr_trm_dgain_s0_ch1, gyr_trm_dgain_s0_ch2, gyr_trm_dgain_s0_ch3] = calculate_sens(obj)
%% CALCULATE_SENS Calculate Gyro Sensitivity
% Gyro Sensitivity  is defined as the slope relating input angular rate  to 
% digital output  in LSB/dps 
% 
% $$\begin{array}{rcl}R & = & S \cdot \Omega \\S & = & S_{MEMS} \cdot S_{ASIC} 
% \cdot S_{PHASE} \\& = & \left(S_{DRIVE} \cdot S_{SENSE}\right) \cdot S_{ASIC} 
% \cdot S_{PHASE}\end{array}$$
%% Gyro Drive Sensitivity
% The Coriolis Force can be expressed as 
% 
% $$\begin{array}{rcl}F_{c} & = & 2\,m_{d}\,x_{d}\,\omega_{d} \cdot \Omega \\S_{DRIVE} 
% & = & \partial F_{c} / \partial \Omega \\& = & 2\,m_{d}\,x_{d}\,\omega_{d}\end{array}$$
% 
% In order to evaluate the sensitivity, we need to evaluate  as part of the 
% AGC loop.
% 
% $$\begin{array}{rcl}C_{s}\left(t\right) & = & C_{0}\sin{\left(\omega_{d} t\right)} 
% \\& = & \left(S_{da} \times x_{da}\right)\sin{\left(\omega_{d} t\right)}\end{array}$$
% 
% The sinusoidal variation of the sensing capacitor  generates a current that 
% is integrated on the  capacitors, and assuming that  is constant:
% 
% $$\begin{array}{rcl}I_{s}\left(t\right) & = & C_{s}\left(t\right)\cdot \frac{\partial 
% V_{GM}}{\partial t} + V_{GM} \cdot \frac{\partial C_{s}\left(t\right)}{\partial 
% t} \\& = & -V_{GM} \cdot \omega_{d} \cdot C_{0} \cos{\left(\omega_{d}t\right)}\end{array}$$
% 
% The output voltage of the Drive CV is then given by
% 
% $$\begin{array}{rcl}V_{out,cvd} & = & \frac{1}{C_{cv,d}}\times \int\limits_{0}^{t}{I_{s}\,d\tau}\\& 
% = & V_{GM}\cdot\frac{C_{0}}{C_{cv,d}}\sin{\left(\omega_{d}t\right)}\end{array}$$
% 
% The demodulation and first-stage integration of the signal from CV
% 
% $$\begin{array}{rcl}I_{r} & = & \frac{V_{out,cv,d}}{R_{sig,d}} \\& = & S_{da} 
% \cdot x_{d} \cdot \frac{V_{GM}}{C_{cv,d}\times R_{sig,d}}\sin{\left(\omega_{d}t\right)} 
% \\& = & I_{r0}\sin{\left(\omega_{d}t\right)}\end{array}$$
%% 
% * The Guard phase has a duration $\left(4/192\right)/f_{d}$
% * The Integration phase has a duration $\big(\left(44-4\right)/192\big)/f_{d}$ 
% * The Guard phase has a duration $\big(\left(48-44\right)/192\big)/f_{d}$ 
%% 
% The charge that is integrated during the four quadrants is:
% 
% $$\begin{array}{rcl}V_{sig,int} & = & Q_{tot,sig}/C_{int,d} \\& = & 4 \cdot 
% Q_{1,sig}/C_{int,d} \\& = & 4/C_{int,d}\times\int\limits_{\frac{4}{192}\cdot\frac{2\pi}{\omega_{d}}}^{\frac{44}{192}\cdot\frac{2\pi}{\omega_{d}}}{I_{r0}\sin{\left(\omega_{d}\tau\right)}d\tau} 
% \\& = & -\big[S_{da}\times \frac{x_{d}}{C_{cv,d}}\times \frac{V_{GM}}{R_{sig,d}}\big]\cdot 
% sinc{\left(\frac{1}{4}\right)}\cdot\sin{\left(\frac{5}{24}\pi\right)}\times 
% \frac{1}{f_{d}\cdot C_{int,d}}\end{array}$$
% 
% Similarly, the reference voltage will generate a current over the resistance  
% which is also integrated over the same four quadrants as follows
% 
% 
% 
% The error signal which is the difference between the two paths t the end of 
% an integration phase should be equal to 0.0
% 
% $$\begin{array}{rcl}V_{ref,int} & = & 4\times\frac{Q_{1,sig}}{C_{int,d}} \\& 
% = & 4\times\frac{5}{24}.\frac{1}{R_{ref,d}}.V_{bg}.\frac{1}{f_{d} \times C_{int,d}}\end{array}$$
% 
% $$\begin{array}{rcl}-V_{sig,int} & = & V_{ref,int} \\x_{d} & = & \frac{4/\pi}{sinc{\left(\frac{1}{4}\right)}\,\cdot\,sinc{\left(\frac{5}{24}\right)}}\cdot 
% \left(\frac{C_{cv,d}}{S_{da}}\right)\times\left(\frac{R_{sig,d}}{R_{ref,d}}\right)\times\left(\frac{V_{bg}}{V_{GM}}\right)\end{array}$$
% 
% Thus, the Gyro Drive Sensitivity is given as:
% 
% $$\begin{array}{rcl}S_{DRIVE} & = & 2\,m_{d}\,\frac{4/\pi}{sinc{\left(\frac{1}{4}\right)}\,\cdot\,sinc{\left(\frac{5}{24}\right)}}\cdot 
% \left(\frac{C_{cv,d}}{S_{da}}\right)\times\left(\frac{R_{sig,d}}{R_{ref,d}}\right)\times\left(\frac{V_{bg}}{V_{GM}}\right)\,\omega_{d}\end{array}$$
    md    = obj.gyro.gyro.stat.mc.md;
    xd    = obj.gyro.gyro.stat.mc.xd;
    Cfbd  = obj.as.gyr.config.inf.drive.cv.gyr_cap_cv_d;
    Sda   = obj.gyro.gyro.stat.mc.Sda;
    Rsigd = obj.as.gyr.config.inf.drive.agc.res_sig;
    Rrefd = obj.as.gyr.config.inf.drive.agc.res_ref;
    Vref  = obj.as.gyr.config.inf.ref.bg.vout;
    Vcm   = obj.as.gyr.config.inf.ref.cm.vout;
    wd    = 2*pi*obj.gyro.gyro.stat.mc.fd;
    sdrive = 2.0*md*(4/pi)/(sinc(1/4)*sinc(5/24))*(Cfbd/Sda)*(Rsigd/Rrefd)*(Vref/Vcm)*wd; % N/dps
%% Gyro Sense Sensitivity
% The sense sensitivity is defined as:
% 
% $$\begin{array}{rcl}C_{s} & = & S_{dd}\times G_{SENSE}\left(\omega_{d}\right) 
% \cdot F_{c}\\S_{SENSE} & = &  \partial C_{s}/\partial F_{c} \\& = & S_{dd}\times 
% G_{SENSE}\left(\omega_{d}\right) \\& = &  S_{dd}\times \frac{1}{m_{s}}\frac{1}{\sqrt{\left(\omega_{s}^2-\omega_{d}^2\right)^2+\left(\omega_{d}\omega_{s}/Q_{s}\right)^2}}\end{array}$$
    ms_C  = obj.gyro.gyro.stat.mc.ms_C;
    ms_S  = obj.gyro.gyro.stat.mc.ms_S;
    ms_Z  = obj.gyro.gyro.stat.mc.ms_Z;
    ws_C = 2*pi*obj.gyro.gyro.stat.mc.fs_C;
    ws_S = 2*pi*obj.gyro.gyro.stat.mc.fs_S;
    ws_Z = 2*pi*obj.gyro.gyro.stat.mc.fs_Z;
    
    Qs_C = obj.gyro.gyro.stat.mc.Qs_C;
    Qs_S = obj.gyro.gyro.stat.mc.Qs_S;
    Qs_Z = obj.gyro.gyro.stat.mc.Qs_Z;
    Cfbs_C = obj.as.gyr.config.inf.sense.cv.C.gyr_cap_cv_s;
    Cfbs_S = obj.as.gyr.config.inf.sense.cv.S.gyr_cap_cv_s;
    Cfbs_Z = obj.as.gyr.config.inf.sense.cv.Z.gyr_cap_cv_s;
    Rrefs = obj.as.gyr.config.inf.sense.rate.common.res_fb;
    Rsigs = obj.as.gyr.config.inf.sense.rate.common.res_sig;
    Sdd_C = 2*md*xd*wd/ms_C/(sqrt(((ws_C/wd)^2-1)^2+(ws_C/wd)^2/Qs_C^2)); % in m/dps
    Sdd_S = 2*md*xd*wd/ms_S/(sqrt(((ws_S/wd)^2-1)^2+(ws_S/wd)^2/Qs_S^2)); % in m/dps
    Sdd_Z = 2*md*xd*wd/ms_Z/(sqrt(((ws_Z/wd)^2-1)^2+(ws_Z/wd)^2/Qs_Z^2)); % in m/dps
    Kin  = obj.as.gyr.config.inf.sense.rate.coeff.Kin;
    Kref = obj.as.gyr.config.inf.sense.rate.coeff.Kref;
    ssense_C = obj.gyro.gyro.stat.mc.Sdd_C_Fdps / (2*md*xd*wd); % F/N
    ssense_S = obj.gyro.gyro.stat.mc.Sdd_S_Fdps / (2*md*xd*wd); % F/N
    ssense_Z = obj.gyro.gyro.stat.mc.Sdd_Z_Fdps / (2*md*xd*wd); % F/N
    srate_C = Kin/Kref*(Vcm/Vref)*(1.0/Cfbs_C)*(Rrefs/(Rsigs)); % 1/F
    srate_S = Kin/Kref*(Vcm/Vref)*(1.0/Cfbs_S)*(Rrefs/(Rsigs)); % 1/F
    srate_Z = Kin/Kref*(Vcm/Vref)*(1.0/Cfbs_Z)*(Rrefs/(Rsigs)); % 1/F
    %5658
    cic_gain_corr = obj.as.gyr.calc_cic_gain_corr(4,1,88);
    afe_gain_corr = fi(5874/8000-1.0,1,10,9);
    gyr_trm_dgain_s0_ch1 = fi(((4000/2^23)-(sdrive*ssense_C*srate_C*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*2^24)^(-1))/(4000/2^23),1,11,12);
    gyr_trm_dgain_s0_ch2 = fi(((4000/2^23)-(sdrive*ssense_S*srate_S*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*2^24)^(-1))/(4000/2^23),1,11,12);
    gyr_trm_dgain_s0_ch3 = fi(((4000/2^23)-(sdrive*ssense_Z*srate_Z*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*2^24)^(-1))/(4000/2^23),1,11,12);
    stot_C = (sdrive*ssense_C*srate_C*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*(1-double(gyr_trm_dgain_s0_ch1))*2^24)^(-1); 
    stot_S = (sdrive*ssense_S*srate_S*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*(1-double(gyr_trm_dgain_s0_ch2))*2^24)^(-1); 
    stot_Z = (sdrive*ssense_Z*srate_Z*(1.0+double(afe_gain_corr))*(1.0+cic_gain_corr.float)*(1-double(gyr_trm_dgain_s0_ch3))*2^24)^(-1); 
end
