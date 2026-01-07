%    so.fsmin  = 8*64*so.se_fd;
%    so.reltol = 1e-9;
%    so.abstol = 1e-11;
function init_gyr_noise(obj)
    %%  noise sampling time and duration; temporary variables
    obj.noise.fsampling = obj.config.inf.drive.vco.gyro_drv_fres*64;
    obj.noise.en = 0;
    
    %%  enabling of noise sources in supply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % common                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bias.VREF_PSD                       = 700e-9^2; %[V^2/Hz]
    %bias.VREF_noise_BW                  = 20e3;     %[Hz]
    %bias.VREF_PSD                       = 257e-9^2; %[V^2/Hz]
    %bias.VREF_noise_BW                  = 75e3;     %[Hz]
    %bias.VREF_PSD_1                     = (12e-6*2)^2; %[V^2/Hz]
    %bias.VREF_noise_BW_1                = 10e-3;     %[Hz]
    %bias.CM_PSD                         = 2e-6^2;   %[V^2/Hz]
    %bias.CM_noise_BW                    = 10e3;     %[Hz]
    
    obj.cm.noise.vcm.en              = and(1, obj.noise.en);
    obj.cm.noise.vref.en             = and(1, obj.noise.en);
    obj.cm.noise.vref.PSD            = 257e-9^2; %[V^2/Hz]
    obj.cm.noise.vref.bw             = 75e3;     %[Hz]
    obj.cm.noise.vcm.PSD             = 2e-6^2;   %[V^2/Hz]
    obj.cm.noise.vcm.bw              = 10e3;     %[Hz]
    
    %%  enabling of noise sources in sense channel
    obj.sense.fe.cv.noise.en                 = and(1, obj.noise.en);
    obj.sense.fe.rate.noise.ampl.en          = and(1, obj.noise.en);
    obj.sense.fe.rate.noise.resin.en         = and(1, obj.noise.en);
    obj.sense.fe.rate.noise.resfb.en         = and(1, obj.noise.en);
    obj.sense.fe.rate.noise.resqc.en         = and(1, obj.noise.en);
    obj.sense.fe.rate.noise.resquadint.en    = and(1, obj.noise.en);
    
    obj.sense.fe.quad.noise.ampl.en          = and(1, obj.noise.en);
    obj.sense.fe.quad.noise.resin.en         = and(1, obj.noise.en);
    obj.sense.fe.quad.noise.resfb.en         = and(1, obj.noise.en);
    
    %%  enabling of noise sources in drive channel
    obj.drive.fe.cv.noise.en         = and(1, obj.noise.en);
    
    obj.drive.fe.agc.noise.ampl.en   = and(1, obj.noise.en);
    obj.drive.fe.agc.noise.resin.en  = and(1, obj.noise.en);
    obj.drive.fe.agc.noise.resfb.en  = and(1, obj.noise.en);
    
    obj.drive.noise.vaa.en           = and(1, obj.noise.en);
    obj.drive.noise.vai.en           = and(1, obj.noise.en);
    
    obj.drive.fe.pll.noise.ampl.en   = and(1, obj.noise.en);
    obj.drive.fe.pll.noise.resin.en  = and(1, obj.noise.en);
    obj.drive.fe.pll.noise.resfb.en  = and(1, obj.noise.en);
    
    obj.drive.be.pll.noise.out.en    = and(1, obj.noise.en);
    obj.drive.be.pll.noise.in.en     = and(1, obj.noise.en);
    
    obj.drive.vco.noise.ressig.en    = and(1, obj.noise.en);
    obj.drive.vco.noise.resref.en    = and(1, obj.noise.en);
    obj.drive.vco.noise.pout.en      = and(1, obj.noise.en);
    
    %%  noise in sense channel
    % charge amplifier noise x/y/z-channels
    %obj.sense.fe.cv.noise.ASD  = [13.7e-9 13.7e-9 13.7e-9]; % amplifier noise ASD [V/sqrt(Hz)]
    %obj.sense.fe.cv.noise.ASD  = [12.4e-9 12.4e-9 12.4e-9]; % amplifier noise ASD [V/sqrt(Hz)]
    obj.sense.fe.cv.noise.ASD  = 13.0e-9; % amplifier noise ASD [V/sqrt(Hz)]
    
    obj.sense.fe.cv.noise.fS   = obj.noise.fsampling;
    obj.sense.fe.cv.noise.seed = 29124;
    
    % integrator noise in rate x/y/z-channels
    % amplifier noise (the same for all channels)
    obj.sense.fe.rate.noise.ampl.ASD  = 35e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.sense.fe.rate.noise.ampl.fS   = obj.noise.fsampling;
    obj.sense.fe.rate.noise.ampl.seed = 23568;
    
    % resistors noise (the same for all channels)
    obj.sense.fe.rate.noise.resin.fS   = obj.noise.fsampling;
    obj.sense.fe.rate.noise.resin.seed = 27658;
    
    obj.sense.fe.rate.noise.resfb.fS   = obj.noise.fsampling;
    obj.sense.fe.rate.noise.resfb.seed = 25854;
    
    obj.sense.fe.rate.noise.resqc.fS   = obj.noise.fsampling;
    obj.sense.fe.rate.noise.resqc.seed = 28653;
    
    obj.sense.fe.rate.noise.resquadint.fS   = obj.noise.fsampling;
    obj.sense.fe.rate.noise.resquadint.seed = 22657;
    
    % integrator noise in quad x/y/z-channels
    obj.sense.fe.quad.noise.ampl.ASD  = 45e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.sense.fe.quad.noise.ampl.fS   = obj.noise.fsampling;
    obj.sense.fe.quad.noise.ampl.seed = 24568;
    
    % resistors noise (the same for all channels)
    obj.sense.fe.quad.noise.resin.fS   = obj.noise.fsampling;
    obj.sense.fe.quad.noise.resin.seed = 26988;
    
    obj.sense.fe.quad.noise.resfb.fS   = obj.noise.fsampling;
    obj.sense.fe.quad.noise.resfb.seed = 25465;
    
    obj.sense.fe.quad.noise.resqc.fS   = obj.noise.fsampling;
    obj.sense.fe.quad.noise.resqc.seed = 28658;
    
    obj.sense.fe.quad.noise.resquadint.fS   = obj.noise.fsampling;
    obj.sense.fe.quad.noise.resquadint.seed = 22659;
    
    %%  noise in drive channel
    % charge amplifier noise
    obj.drive.fe.cv.noise.ASD  = obj.config.inf.drive.cv.noise.Vn_opa;
    obj.drive.fe.cv.noise.fS   = obj.noise.fsampling;
    obj.drive.fe.cv.noise.seed = 82457;
    % resistor noise
    
    %% integrator noise in agc loop
    % amplifier noise
    obj.drive.fe.agc.noise.ampl.ASD  = 35e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.drive.fe.agc.noise.ampl.fS   = obj.noise.fsampling;
    obj.drive.fe.agc.noise.ampl.seed = 80134;
    % resistors noise
    obj.drive.fe.agc.noise.resin.fS   = obj.noise.fsampling;
    obj.drive.fe.agc.noise.resin.seed = 87658;
    
    obj.drive.fe.agc.noise.resfb.fS   = obj.noise.fsampling;
    obj.drive.fe.agc.noise.resfb.seed = 83659;
    
    %% integrator noise in pll loop
    % amplifier noise
    obj.drive.fe.pll.noise.ampl.ASD  = 38e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.drive.fe.pll.noise.ampl.fS   = obj.noise.fsampling;
    obj.drive.fe.pll.noise.ampl.seed = 90134;
    
    % resistors noise
    obj.drive.fe.pll.noise.resin.fS   = obj.noise.fsampling;
    obj.drive.fe.pll.noise.resin.seed = 97658;
    
    obj.drive.fe.pll.noise.resfb.fS   = obj.noise.fsampling;
    obj.drive.fe.pll.noise.resfb.seed = 93659;
    
    % be input noise
    obj.drive.be.pll.noise.in.ASD  = 500e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.drive.be.pll.noise.in.fS   = 70e3;
    obj.drive.be.pll.noise.in.seed = 94136;
    
    % be output noise
    obj.drive.be.pll.noise.out.ASD  = 2000e-9; % amplifier noise ASD [V/sqrt(Hz)]
    obj.drive.be.pll.noise.out.fS   = 140e3;
    obj.drive.be.pll.noise.out.seed = 95673;
    
    
    % vco noise
    % resistors noise
    obj.drive.vco.noise.ressig.fS   = obj.noise.fsampling;
    obj.drive.vco.noise.ressig.seed = 97623;
    
    obj.drive.vco.noise.resref.fS   = obj.noise.fsampling;
    obj.drive.vco.noise.resref.seed = 91689;
    
    obj.drive.vco.noise.pout.ASD    = 2000e-9; % pout noise ASD [V/sqrt(Hz)]
    obj.drive.vco.noise.pout.fS     = obj.noise.fsampling;
    obj.drive.vco.noise.pout.seed   = 92304;
end
