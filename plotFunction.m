function data = plotFunction(varargin)

% Check & assigning arguments
% Number of input arguments
nin = nargin;
% Number of output arguments
nargoutchk(0,5)

% Checking arguments
data = ChkArgs(varargin,nin);

if ~data.flags.forcedirectplot % Draw plots
    data = FncSignalSizeCorrection(data);
    if data.flags.SpectralCrossCorrelation 
        if data.flags.calculate_IQ
            data = FncPredictPhiDemodNoise(data);
        end
        switch 2 % Select function to predict spectral correlation
            case 1
                data = FncSpectralCrossCorrelation_IBN_Fast(data);
            case 2
                data = FncSpectralCrossCorrelation_Rate_Slow(data);
            case 3
                data = FncSpectralCrossCorrelation_Rate_Fast(data);
        end
    elseif data.flags.DigitalCorrection
        data = FncPowerSpecIQ(data);
        data = FncDigitalCorrection(data);
    else
        
        % Calculate spectrum
        data = FncPowerSpec(data);
        
        % Predict Bins
        data = FncCalculateBINs(data);
        
        % Calculate signal power
        data = FncSignalPower(data);
        
        % Calculate IBN
        data = FncEstimatingIBN(data);
        
        % calculate avarage noise density with defined frequency range
        data = FncPredictAvarageNoiseDensity(data);
        
        if data.flags.calculate_IQ
            data = FncPowerSpecIQ(data);
        end
    end
else
    data.flags.new_fig = true;
end

% Call plot function
if data.flags.plot_fft
    data = FncPlotPSD(data);
end

if data.info.verbose
    fprintf(data.info.txt)
end

if data.flags.ReducedMemoryStorage
    data = FncReduceMemoryStorage(data);
end

function data = FncReduceMemoryStorage(data)
data.timesignal.v = [];
data.timesignal.v_I = [];
data.timesignal.v_Q = [];
data.spectrum.W = [];

function data = FncSignalSizeCorrection(data)
% SignalSizeCorrection -------------------- ['ForceTruncationPower2'|'ForceTruncationFactor2'|'ForcePaddingPower2','ForcePaddingFactor2']
% SignalSizeCorrectionMethod -------------- help padarray ['zeros','circular','replicate','symmetric']
v = data.timesignal.v;
Np = size(v,1);
% Check if vector length is power of two

switch data.timesignal.parameters.SignalSizeCorrection
    case {'ForceTruncationPower2'}
        if mod(log2(Np),1) % This is more strict compared to the requirements of the fft
            Npm = 2^(nextpow2(Np)-1);
            v = v(1:Npm,1);
            fprintf('plotFunction:input signal was truncated.\n')
        end
    case {'ForceTruncationFactor2'}
        if mod(Np,2) % This is more strict compared to the requirements of the fft
            Npm = Np-1;
            v = v(1:Npm,1);
            fprintf('plotFunction:input signal was truncated.\n')
        end
    case {'ForcePaddingPower2'}
        if mod(log2(Np),1) % This is more strict compared to the requirements of the fft
            Npm = 2^(nextpow2(Np));
            v = v(1:Npm,1);
            if strcmpi(data.timesignal.parameters.SignalSizeCorrectionMethod,'zeros')
                v = padarray(v,[Npm 0],0,'post');
            else
                v = padarray(v,[Npm 0],data.timesignal.parameters.SignalSizeCorrectionMethod,'post');
            end
            fprintf('plotFunction:input signal was padded.\n')
        end
    case {'ForcePaddingFactor2'}
        if mod(log2(Np),1) % This is more strict compared to the requirements of the fft
            Npm = Np+1;
            if strcmpi(data.timesignal.parameters.SignalSizeCorrectionMethod,'zeros')
                v = padarray(v,[Npm 0],0,'post');
            else
                v = padarray(v,[Npm 0],data.timesignal.parameters.SignalSizeCorrectionMethod,'post');
            end
            fprintf('plotFunction:input signal was padded.\n')
        end
    otherwise
        error('plotFunction:SignalSizeCorrection','Defined option is not supported.')
end

data.timesignal.v = v;
data.timesignal.parameters.Np = size(v,1);
data.spectrum.parameters.segment_length = floor(data.timesignal.parameters.Np/data.spectrum.parameters.segment_no);
data.spectrum.parameters.overlap_no = floor(data.timesignal.parameters.Np/data.spectrum.parameters.segment_no*data.spectrum.parameters.overlapRatio);
data.spectrum.parameters.overlap_string = floor(data.spectrum.parameters.overlapRatio * 100);

function data = FncSpectralCrossCorrelation_IBN_Fast(data)
v = data.timesignal.v;
NOVERLAP = data.spectrum.parameters.overlap_no;
NOPOINTS = data.spectrum.parameters.segment_length;
no_nz_bins = data.spectrum.parameters.no_nz_bins;
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;

% Segmenting the input signal
nx = length(v);
ncol = fix((nx-NOVERLAP)/(NOPOINTS-NOVERLAP));
colindex = 1 + (0:(ncol-1))*(NOPOINTS-NOVERLAP);
rowindex = (1:NOPOINTS)';
vin = zeros(NOPOINTS,ncol);
vin(:) = v(rowindex(:,ones(1,ncol))+colindex(ones(NOPOINTS,1),:)-1);

W = blackman(NOPOINTS,'periodic');

Pxy_2s_a = zeros(NOPOINTS,ncol);

if ~mod(NOPOINTS,2)
    nsel = NOPOINTS/2;
else
    nsel = ceil(NOPOINTS/2);
end

switch data.options.CorrelationFunction
    case 'cxcorr'
        hCF = @circconv;
    case 'xcorr'
        hCF = @classconv;
    otherwise
        error('plotFunction:CrossCorrelation','Correlation method is not supported');
end

f = (0:nsel-1)/NOPOINTS*fs;
sel = f>fs/2;
f(sel) = f(sel)-fs;
posf0 = find(f>f0,1,'first');

nfb = data.options.FocusInBand;
if nfb
    fbsel = (f>f0-nfb*fb)&(f<f0+nfb*fb);
    if posf0 < 3
        fbsel(1:3) = false;
    else
        fbsel(posf0-2:posf0+2) = false;
    end
    parfor ni = 1:ncol
        Fv = (fft(vin(:,ni).*W))/sum(W);
        Fv2 = Fv;
        Fv2(~fbsel) = 0;
        c = xcov(Fv,Fv2,'coeff');
        Pxy_2s_a(:,ni) = abs(c(1:NOPOINTS));
    end
else
    for ni = 1:ncol
        W = 1;
        Fv = (fft(vin(:,ni).*W))/sum(W);
        Fv(posf0-no_nz_bins:posf0+no_nz_bins)=0;  % Cancle signals at positive frequency
        Fv(end-2-(posf0+no_nz_bins):end-2-(posf0-no_nz_bins))=0; % Cancle signals at negative frequency
        c = hCF(Fv);
        Pxy_2s_a(:,ni) = abs(c(1:NOPOINTS)).^2;
    end
end
% Pxy_re = abs(real(c(1:dimv/2)));
% Pxy_im = abs(imag(c(1:dimv/2)));
% Pxy_abs = abs(c(1:dimv/2));
Pxy_2s = sum(Pxy_2s_a,2);

Pxy = Pxy_2s(1:nsel);
Pxy(1:5) = NaN;

Np = NOPOINTS;
data.spectrum.psd = sqrt(Pxy);
data.timesignal.parameters.Np = Np; % Segment length changes the number of points!
data.spectrum.f = f;
data.spectrum.W = W;
data.spectrum.NG = sum(W.^2)/Np;
data.spectrum.CG = sum(W)/Np;
data.spectrum.dfbin = fs/Np;
% Set structure for plot
data.spectrum.psd_IB_accumulated = [];
data.spectrum.parameters.bins_fband = [];
data.spectrum.parameters.bins_signal = [];

function data = FncSpectralCrossCorrelation_Rate_Slow(data)
V = data.timesignal.v;
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;
FixPhase = data.timesignal.parameters.phi_demod;
FS = data.spectrum.parameters.FS;
no_nz_bins = data.spectrum.parameters.no_nz_bins;
FreqMax = data.options.MaximumFrequency;
dFreq = data.options.FrequencyResolution;

Np = size(V,1);

W = blackman(Np,'periodic'); % Window data to reduce noise leakage
vW = FS/2*V.*W;
Wnorm = sqrt(W'*W*Np);% norm(W,2)*sqrt(dimv/2);
Vc = fft(vW); % TODO: Norm is not right. Has to be corrected. Done.
Vc(2:end-1) = Vc(2:end-1)*sqrt(2);
Vc = Vc/Wnorm;
Vc = [Vc; Vc];
dfnorm = fs/Np;
% fnorm = (-Np:Np-1)*dfnorm;
Nfbnorm = round(fb/fs*Np);
Nf0norm = round(f0/fs*Np);

% fb = round(fb/fs*Np)*fs/Np;

Vref = FncModSpectrumN(Vc,Np,0,Nf0norm,Nfbnorm,pi/2,-FixPhase);
Vref(1:no_nz_bins)=[];

if size(Vref,1)<5
    errordlg('Not enought simulation points for correlation','Spectral correlation')
    return
end

r = NaN(Np/2+1,1);

if Nf0norm > 0
    fpoints = Nf0norm*(1:1:5);
else
    fpoints = Np/16*(1:1:8);
end

if FreqMax > fs/2
    FreqMax = fs/2;
end

fnp = union(0:round(dFreq/dfnorm):round(FreqMax/dfnorm),fpoints);
% check if you will be above fs/2
fnp(fnp>Np/2)=[];

Nc = size(fnp,2);
rtI = NaN(Nc,1);
rtQ = NaN(Nc,1);
% pvallim = 0.10;
for nc = 1:Nc
    ni = fnp(nc);
    Vdem = FncModSpectrumN(Vc,Np,ni,Nf0norm,Nfbnorm,-FixPhase,-FixPhase);
    Vdem(1:no_nz_bins)=[];
    tmp = corr(Vref,Vdem);
%     [tmp, pval] = corr([real(Vref);imag(Vref)],[real(Vdem);imag(Vdem)]);
%     if pval > pvallim
%         tmp = NaN;
%     end
    rtI(nc) = norm(tmp);
    Vdem = FncModSpectrumN(Vc,Np,ni,Nf0norm,Nfbnorm,-FixPhase+pi/2,-FixPhase);
    Vdem(1:no_nz_bins)=[];
    tmp = corr(Vref,Vdem);
%     [tmp, pval] = corr([real(Vref);imag(Vref)],[real(Vdem);imag(Vdem)]);
%     if pval > pvallim
%         tmp = NaN;
%     end
    rtQ(nc) = norm(tmp);
    %     if (ni-tni)/fnpe >= 0.1
    %         fprintf('.')
    %         tni = ni;
    %     end
end
% fprintf('Done.\n')
r(fnp+1) = rtI+1i*rtQ;

% INFO: Can this calculation be done by using cconv? cconv function is much
% faster.

data.spectrum.psd = r;
% data.timesignal.parameters.Np = Np; % Segment length changes the number of points!
data.spectrum.f = (0:Np/2)*fs/Np;
data.spectrum.W = W;
data.spectrum.NG = sum(W.^2)/Np;
data.spectrum.CG = sum(W)/Np;
data.spectrum.dfbin = fs/Np;
% Set structure for plot
data.spectrum.psd_IB_accumulated = [];
data.spectrum.parameters.bins_fband = [];
data.spectrum.parameters.bins_signal = [];

function data = FncSpectralCrossCorrelation_Rate_Fast(data)
V = data.timesignal.v;
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;
FixPhase = data.timesignal.parameters.phi_demod;
FS = data.spectrum.parameters.FS;
no_nz_bins = data.spectrum.parameters.no_nz_bins;

Np = size(V,1);
W = blackman(Np,'periodic');
vW = FS/2*V.*W;
Wnorm = sqrt(W'*W*Np);% norm(W,2)*sqrt(dimv/2);
Vc = fft(vW); % TODO: Norm is not right. Has to be corrected.
Vc(2:end-1) = Vc(2:end-1)*sqrt(2);
Vc = Vc/Wnorm;
f = (0:Np-1)*fs/Np;
sel = f>fs/2;
f(sel) = f(sel)-fs;
df = fs/Np;

Vref = FncModSpectrum(Vc,f,0,f0,fb,pi/2,-FixPhase);
Vref(1:no_nz_bins)=[];

phase_rad_pp = -FixPhase + pi/2;
phase_rad_pn = -FixPhase - pi/2;

vph_pp = exp(1i*phase_rad_pp);
vph_nn = conj(vph_pp);
vph_pn = exp(1i*phase_rad_pn);
vph_np = conj(vph_pn);

fsel = -f0+f;
sel = find(fsel==0);
Vc_pp = circshift(Vc,sel-1);

fsel = -f0+f;
sel = find(fsel==0);
Vc_nn = circshift(Vc,sel-1);

fsel = -f0+f;
sel = find(fsel==0);
Vc_pp = circshift(Vc,sel-1);

Vcr_pp = 1i*vph_nn*Vc_pp;
Vcr_nn = -1i*vph_pp*Vc_nn;
Vcr_pn = -1i*vph_np*Vc_pn;
Vcr_np = 1i*vph_pn*Vc_np;

data.spectrum.psd = r;
data.timesignal.parameters.Np = Np; % Segment length changes the number of points!
data.spectrum.f = f(~sel);
data.spectrum.W = W;
data.spectrum.NG = sum(W.^2)/Np;
data.spectrum.CG = sum(W)/Np;
data.spectrum.dfbin = fs/Np;
% Set structure for plot
data.spectrum.psd_IB_accumulated = [];
data.spectrum.parameters.bins_fband = [];
data.spectrum.parameters.bins_signal = [];


function c = circconv(u)
c = ifft(fft(u).*conj(fft(u)));

function c = classconv(u)
c = xcov(u,u,'coeff');

function data = FncPowerSpec(data)
% Extract parameters
v = data.timesignal.v;
segment_length = data.spectrum.parameters.segment_length;
overlap_no = data.spectrum.parameters.overlap_no;
fs = data.timesignal.parameters.fs;

% Check for NaN or invalid data
if isempty(v) || all(isnan(v))
    error('FncPowerSpec:InvalidInput', 'Input signal v contains only NaN values or is empty. Check upstream simulation.');
end

% Remove NaN values if they exist
if any(isnan(v))
    warning('FncPowerSpec:NaNDetected', 'Input signal contains NaN values. Removing NaN entries before pwelch.');
    v = v(~isnan(v));
    if length(v) < segment_length
        error('FncPowerSpec:InsufficientData', 'After removing NaN values, insufficient data remains (need at least %d samples, have %d).', segment_length, length(v));
    end
end

% Windowing
W=blackman(segment_length,'periodic');

[ys,f] = pwelch(v,W,overlap_no,segment_length,fs,'onesided');

Np = 2*length(ys); % Segment length changes the number of points! Recalculate effective number of points assuming that the FFT halfs the number of points.
data.spectrum.psd = ys;
data.timesignal.parameters.Np = Np; % Segment length changes the number of points!
data.spectrum.f = f;
data.spectrum.W = W;
data.spectrum.NG = sum(W.^2)/Np;
data.spectrum.CG = sum(W)/Np;
data.spectrum.dfbin = fs/Np;

function data = FncPowerSpecIQ(data)
% Get data &  parameters
V = data.timesignal.v;
dimv = length(V);
segment_length = data.spectrum.parameters.segment_length;
overlap_no = data.spectrum.parameters.overlap_no;
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;
no_nz_bins = data.spectrum.parameters.no_nz_bins;

if data.flags.flag_predict_phi_demod
    if strcmpi(data.options.predict_phi_demod,'N')
        data = FncPredictPhiDemodNoise(data);
    else
        data = FncPredictPhiDemod(data);
    end
end

phi_demod = data.timesignal.parameters.phi_demod;
dfbin = data.spectrum.dfbin;

% Comment:
% I-channel is real, Q-channel is imaginary
% TODO: Use Frequency components to calculate I and Q spectrum.
v_I = V.*cos(2*pi*f0/fs*(0:dimv-1)-phi_demod)'*2; % real component
v_Q = V.*sin(2*pi*f0/fs*(0:dimv-1)-phi_demod)'*2; % imaginary component
W = data.spectrum.W;

[psd_I,~]= pwelch(v_I,W,overlap_no,segment_length,fs,'onesided');
[psd_Q,f]= pwelch(v_Q,W,overlap_no,segment_length,fs,'onesided');


fpos = find(f<=10*fb,1,'last');
fband = (f<=fb);
switch lower(data.spectrum.parameters.iq_noise)
    case 'two_sided'
        data.spectrum.psd_I = [psd_I(fpos:-1:2)/2; psd_I(1); psd_I(2:fpos)/2]; % Two sided spectrum
        data.spectrum.psd_Q = [psd_Q(fpos:-1:2)/2; psd_Q(1); psd_Q(2:fpos)/2]; % Two sided spectrum
    case 'one_sided'
        data.spectrum.psd_I = [psd_I(fpos:-1:2); psd_I(1); psd_I(2:fpos)]; % One sided spectrum
        data.spectrum.psd_Q = [psd_Q(fpos:-1:2); psd_Q(1); psd_Q(2:fpos)]; % One sided spectrum
    otherwise
        error('plotFunction:iq_noise','%s not supported',lower(data.spectrum.parameters.iq_noise));
end
data.spectrum.f_IQ = f0+[-f(fpos:-1:2); f(1); f(2:fpos)];
P_sig_I = sum(psd_I(1:no_nz_bins))*dfbin;
P_sig_Q = sum(psd_Q(1:no_nz_bins))*dfbin;
data.signal.phi_IQ = atan2(P_sig_Q,P_sig_I);
data.signal.phi_IQ_deg = data.signal.phi_IQ*180/pi;

% Remove signal from spectrum
psd_I(1:no_nz_bins) = psd_I(no_nz_bins+1);
psd_Q(1:no_nz_bins) = psd_Q(no_nz_bins+1);

P_sig_I_norm = P_sig_I/data.spectrum.PowerFS;
P_sig_Q_norm = P_sig_Q/data.spectrum.PowerFS;
P_sig_I_renorm = P_sig_I_norm*data.spectrum.PowerRescaleFS;
P_sig_Q_renorm = P_sig_Q_norm*data.spectrum.PowerRescaleFS;

P_IBN_I = sum(psd_I(fband))*dfbin;
P_IBN_Q = sum(psd_Q(fband))*dfbin;
P_IBN_I_norm = P_IBN_I/data.spectrum.PowerFS;
P_IBN_I_renorm = P_IBN_I_norm*data.spectrum.PowerRescaleFS;
P_IBN_Q_norm = P_IBN_Q/data.spectrum.PowerFS;
P_IBN_Q_renorm = P_IBN_Q_norm*data.spectrum.PowerRescaleFS;

data.signal.P_signal_I = P_sig_I;
data.signal.P_signal_Q = P_sig_Q;
data.signal.P_signal_I_norm = P_sig_I_norm;
data.signal.P_signal_Q_norm = P_sig_Q_norm;
data.signal.P_signal_I_renorm = P_sig_I_renorm;
data.signal.P_signal_Q_renorm = P_sig_Q_renorm;
data.ibn.P_IBN_I = P_IBN_I;
data.ibn.P_IBN_Q = P_IBN_Q;
data.ibn.P_IBN_I_norm = P_IBN_I_norm;
data.ibn.P_IBN_I_renorm = P_IBN_I_renorm;
data.ibn.P_IBN_Q_norm = P_IBN_Q_norm;
data.ibn.P_IBN_Q_renorm = P_IBN_Q_renorm;

data.timesignal.v_I = v_I;
data.timesignal.v_Q = v_Q;

% function [Vt,ncol] = SegmentingBitstream(V,NOVERLAP,NOPOINTS)
% nx = length(V);
% ncol = fix((nx-NOVERLAP)/(NOPOINTS-NOVERLAP));
% colindex = 1 + (0:(ncol-1))*(NOPOINTS-NOVERLAP);
% rowindex = (1:NOPOINTS)';
% Vt = zeros(NOPOINTS,ncol);
% Vt(:) = V(rowindex(:,ones(1,ncol))+colindex(ones(NOPOINTS,1),:)-1);

function data = FncPredictPhiDemod(data)
V = data.timesignal.v;
dimv = length(V);
% dimv = 2^(nextpow2(dimv)-1);
% V = V(1:dimv);
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;

v_I = V.*cos(2*pi*f0/fs*(0:dimv-1))'*2;
v_Q = V.*sin(2*pi*f0/fs*(0:dimv-1))'*2;
phi_demod = atan2(mean(v_Q),mean(v_I));
if strcmpi(data.options.predict_phi_demod,'i')
    data.timesignal.parameters.phi_demod = phi_demod;
else
    data.timesignal.parameters.phi_demod = phi_demod+pi/2;
end

function data = FncPredictPhiDemodNoise(data)
% Under the assumption that the rate noise is much smaller compared to the
% quadrature noise, the phase can also be predicted by searching for the
% minimum of the rate IBN.
V = data.timesignal.v;
dimv = length(V);
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;
no_nz_bins = data.spectrum.parameters.no_nz_bins;

W = blackman(dimv,'periodic');
Wnorm = norm(W,2)*sqrt(dimv/2);
Npbw = round(dimv/fs*fb);
vW = V.*W;
Vc = fft(vW)/Wnorm;
f = (0:dimv-1)*fs/dimv;
sel = f>fs/2;
f(sel) = f(sel)-fs;
pos = find(f>=f0,1,'first');

Vcp = Vc(pos+(1:Npbw));
Vcn = Vc(end-(pos-2-(1:Npbw)));

% Coars optimization to find the minimum rate IBN
parray = pi/180*(-180:1:0);
Npa = size(parray,2);
IBNv = zeros(Npa,1);
% figure;ha = axes('next','add','YScale','log');
parfor ni = 1:Npa
    phase_rad = parray(ni);
    vphp = exp(1i*phase_rad);
    vphn = conj(vphp);
    Vcpp = vphp*Vcp;
    Vcnn = vphn*Vcn;
    Vdem = (Vcnn+Vcpp); % Factor of 1/2 ommited as it cancles with the factor of 2 from the demodulation %     Vdem = ModSpectrum(Vc,f,0,f0,fb,0,phase_rad); % Not used due to speed
    Vdem(1:no_nz_bins)=0;
%     plot(ha,abs(Vdem))
    IBNv(ni) = Vdem'*Vdem;
end
[~,ind]=min(IBNv);

if ind-1 < 1
    ind = ind+1;
elseif ind+1 > Npa
    ind = ind-1;
end

% Second run to get a better estimation
parray = linspace(parray(ind-1),parray(ind+1),data.options.finephaseno+1);
Npa = size(parray,2);
IBNv = zeros(Npa,1);

parfor ni = 1:Npa
    phase_rad = parray(ni);
    vphp = exp(1i*phase_rad);
    vphn = conj(vphp);
    Vcpp = vphp*Vcp;
    Vcnn = vphn*Vcn;
    Vdem = (Vcnn+Vcpp); % Factor of 1/2 ommited as it cancles with the factor of 2 from the demodulation
    Vdem(1:no_nz_bins)=0;
    IBNv(ni) = Vdem'*Vdem;
end
[~,ind]=min(IBNv);

data.timesignal.parameters.phi_demod = parray(ind);
data.spectrum.IQ.Vcp = Vcp;
data.spectrum.IQ.Vcn = Vcn;
data.spectrum.IQ.exp_i_phi =  exp(parray(ind));

function data = FncDigitalCorrection(data)
% Warning: this function is not used and currently not supported
% Under the assumption that the rate noise is much smaller compared to the
% quadrature noise, the phase can also be predicted by searching for the
% minimum of the rate IBN.
V = data.timesignal.v;
fs = data.timesignal.parameters.fs;
f0 = data.timesignal.parameters.f0;
fb = data.timesignal.parameters.fb;
fmodv = data.spectrum.parameters.fmod;
FixPhase = data.timesignal.parameters.phi_demod;
FS = data.spectrum.parameters.FS;
no_nz_bins = data.spectrum.parameters.no_nz_bins;

[Vdem, ~] = FncDemSpectrum(V,FS,fs,f0,fb,FixPhase,fmodv,no_nz_bins);
data.spectrum.Vdem = Vdem;

function [VdemI, VdemQ] = FncDemSpectrum(V,FS,fs,f0,Bandwith,FixPhase,fmodv,no_nz_bins)
dimv = size(V,1); % Number of input samples
W = blackman(dimv,'periodic'); % Calculate window
vW = FS/2*V.*W; % Weighted signal by window
Wnorm = sqrt(W'*W*dimv);% norm(W,2)*sqrt(dimv/2);
Vc = fft(vW); % Calculate fft
Vc(2:end-1) = Vc(2:end-1)*sqrt(2);
Vc = Vc/Wnorm; % Normalize spectrum
% Vc = fft(vW)/sqrt(sum(W.^2)*fs);
f = (0:dimv-1)*fs/dimv; % Generate frequency vector 0..fs
sel = f>fs/2; % slect the frequency which are above fs/2
f(sel) = f(sel)-fs; % shift 0..fs => -fs/2..fs/2

Nfb = ceil(Bandwith/fs*dimv); % Number of samples within bandwidth

dim_fmod = size(fmodv,2); % Allocated memory for used 
VdemI = zeros(Nfb,dim_fmod); % Allocated memory of demodulated spectrum: rate proportional
VdemQ = zeros(Nfb,dim_fmod); % Allocated memory of demodulated spectrum: quadrature proportional

ni = 1; % Calculate normal spectrum
VdemI(:,ni) = FncModSpectrum(Vc,fnorm,fmodv(ni)*f0norm,f0norm,fbnorm,pi/2,-FixPhase); % TODO: Check arguments. Phase seems to be wrong! Done. Now it is o.k.
VdemI(1:no_nz_bins,ni)=0; % Delete samples close to DC due to leakage of signal power
VdemQ(:,ni) = FncModSpectrum(Vc,fnorm,fmodv(ni)*f0norm,f0norm,fbnorm,0,-FixPhase);
VdemQ(1:no_nz_bins,ni)=0; % Delete samples close to DC due to leakage of signal power

for ni = 2:dim_fmod
    VdemI(:,ni) = FncModSpectrum(Vc,f,fmodv(ni)*f0,f0,Bandwith,-FixPhase,-FixPhase);
    VdemI(1:no_nz_bins,ni)=0;
    VdemQ(:,ni) = FncModSpectrum(Vc,f,fmodv(ni)*f0,f0,Bandwith,-FixPhase,-FixPhase+pi/2);
    VdemQ(1:no_nz_bins,ni)=0;
end

function Vdem = FncModSpectrumN(Vc,Np,Nfmod,Nfdem,Nfb,phi_rad,phi_dem_rad)
Nf1 = Nfdem + Nfmod;
Nf2 = Nfdem - Nfmod;

Vc_pp = Vc(Np+Nf1:Np+Nf1+Nfb); % pp <=> positive fdem; positive fmod
Vc_nn = Vc(Np-Nf1:Np-Nf1+Nfb); % nn <=> negative fdem; negative fmod
Vc_pn = Vc(Np+Nf2:Np+Nf2+Nfb); % pn <=> positive fdem; negative fmod
Vc_np = Vc(Np-Nf2:Np-Nf2+Nfb); % np <=> negative fdem; positive fmod

phase_rad_pp = phi_dem_rad + phi_rad;
phase_rad_pn = phi_dem_rad - phi_rad;

vph_pp = exp(1i*phase_rad_pp);
vph_nn = conj(vph_pp);
vph_pn = exp(1i*phase_rad_pn);
vph_np = conj(vph_pn);

Vcr_pp = 1i*vph_nn*Vc_pp;
Vcr_nn = -1i*vph_pp*Vc_nn;
Vcr_pn = -1i*vph_np*Vc_pn;
Vcr_np = 1i*vph_pn*Vc_np;

% \\rt02fs07.de.bosch.com\MA$\ese\buh2sh\work\Mathematica\20111024_Bm_NoiseModulation.nb
% Vdem = (Vcr_pp(Nminsel,1)+Vcr_nn(Nminsel,1)+Vcr_pn(Nminsel,1)+Vcr_np(Nminsel,1))/2;
Vdem = (Vcr_pp+Vcr_nn+Vcr_pn+Vcr_np)/2;

function data = FncSignalPower(data)
% Extract parameters
psd = data.spectrum.psd;
bins_signal = data.spectrum.parameters.bins_signal;
Nsig = data.spectrum.parameters.Nsig;
dfbin = data.spectrum.dfbin;
PowerFS = data.spectrum.PowerFS;
PowerRescaleFS = data.spectrum.PowerRescaleFS;

% Power of the main signal bin
P_signal = zeros(1,Nsig);
for ni = 1:Nsig
    P_signal(ni) = sum(psd(bins_signal{ni}))*dfbin;
end
P_signal_norm = P_signal / PowerFS;
P_signal_renorm = P_signal_norm * PowerRescaleFS;

data.signal.P_signal = P_signal;
data.signal.P_signal_norm = P_signal_norm;
data.signal.P_signal_renorm = P_signal_renorm;

% Calculation of the signal power in dB
data.signal.P_signal_dBFS = 10*log10(P_signal_norm);
data.signal.P_signal_dB = 10*log10(P_signal);


function data = FncEstimatingIBN(data)
% Extract parameters
psd_noise = data.spectrum.psd;
bins_fband = data.spectrum.parameters.bins_fband;
PowerFS = data.spectrum.PowerFS;
PowerRescaleFS = data.spectrum.PowerRescaleFS;
dfbin = data.spectrum.dfbin;
% bins_fband_intersect = data.spectrum.parameters.bins_fband_intersect;

% TODO: Noise power within no_nz_bins is neglected. This leads to an error,
% if the points wihtin IBN are very small.

% Calculate accumulated noise for each bin
% psd_noise(bins_fband_intersect) = (psd_noise(bins_fband_intersect(1))+psd_noise(bins_fband_intersect(end)))/2;

psd_IB = psd_noise(bins_fband);
% psd_IB_intersect = psd_noise(bins_fband_intersect);
% if length(psd_IB) < 15
%     fprintf('The number of points within the IBN that is used for the calculation is very small.\nThis leads to an error of the predicted IBN\nTo avoide this, increase the number of points wihtin the IBN!\n')
% end
psd_IB_accumulated = cumsum(psd_IB)*dfbin;
if ~isempty(psd_IB_accumulated)
    P_IBN = psd_IB_accumulated(end);
else
    P_IBN = NaN;
end
P_IBN_norm = P_IBN/PowerFS;
P_IBN_renorm = P_IBN_norm * PowerRescaleFS;
data.spectrum.psd_IB = psd_IB;
data.spectrum.psd_IB_accumulated = psd_IB_accumulated;
data.ibn.P_IBN = P_IBN;
data.ibn.P_IBN_norm = P_IBN_norm;
data.ibn.P_IBN_renorm = P_IBN_renorm;
data.ibn.P_IBN_dBFS = 10*log10(P_IBN_norm);
data.ibn.P_IBN_dB = 10*log10(P_IBN);
data.signal.SNDR_dB = data.signal.P_signal_dBFS-data.ibn.P_IBN_dBFS;

function data = FncPredictAvarageNoiseDensity(data)
psd = data.spectrum.psd;
f = data.spectrum.f;
PowerFS = data.spectrum.PowerFS;
PowerRescaleFS = data.spectrum.PowerRescaleFS;
dfbin = data.spectrum.dfbin;
Wnorm = 1;

% Predict noise density within 1...20Hz range
% Predict noise density within 1...100Hz range
fmin = 1;
fmax = 100;
data.spectrum.AverageNoiseDensity.fmin = fmin;
data.spectrum.AverageNoiseDensity.fmax = fmax;
fsel = fmin<=f&f<=fmax;
frange = f(fsel);
% Check that no signal bin can leak into the noise calculation
if find(fsel,1,'first')<=3
    fsel(1:3)=0;
    fprintf('Frequency range for noise density calculation changed\n')
end

if sum(fsel)>5
    Smean = sqrt(sum((psd(fsel)*dfbin/PowerFS*PowerRescaleFS*Wnorm))/(frange(end)-frange(1)));
else
    Smean = NaN;
end
data.spectrum.AverageNoiseDensity.AND = Smean;

function [scale, str] = FncScaleSignal(value)
log_value = log10(value);
if log_value < -6
    scale = 1e9;
    str = 'n';
elseif  log_value < -3
    scale = 1e6;
    str = 'µ';
elseif  log_value < 0
    scale = 1e3;
    str = 'm';
elseif  log_value < 3
    scale = 1;
    str = '';
elseif  log_value < 6
    scale = 1 * 1e-3;
    str = 'k';
elseif  log_value < 9
    scale = 1 * 1e-6;
    str = 'M';
else
    scale = 1 * 1e-9;
    str = 'G';
end

function data = FncPlotPSD(data)
psd = data.spectrum.psd;
psd_IB_accumulated = data.spectrum.psd_IB_accumulated;
new_fig = data.flags.new_fig;
normfs = data.flags.normfs;
scale = data.options.scale;
bins_fband = data.spectrum.parameters.bins_fband;
bins_signal = data.spectrum.parameters.bins_signal;
str_title = data.options.str_title;
flabel = data.flags.flabel;
f = data.spectrum.f;
fs = data.timesignal.parameters.fs;
Np = data.timesignal.parameters.Np;
stats = data.flags.stats;
% W = data.spectrum.W;
NG = data.spectrum.NG;
CG = data.spectrum.CG;
c_signal = data.color.signal;
c_psd = data.color.psd;
c_psd_I = data.color.psd_I;
c_psd_Q = data.color.psd_Q;
c_ibn = data.color.ibn;
PowerFS = data.spectrum.PowerFS;
PowerRescaleFS = data.spectrum.PowerRescaleFS;
dfbin = data.spectrum.dfbin;
show_IQ = data.flags.show_IQ;
fsig = data.timesignal.parameters.fsig;
unit_str = data.timesignal.parameters.unit;
Nsig = data.spectrum.parameters.Nsig;
% New figure
if new_fig
    fprintf('Creating new figure\n')
    hfig = figure;
    clf;
    haxes = axes('Parent',hfig);
elseif ~isa(data.options.haxes,'matlab.graphics.axis.Axes') && isnan(data.options.haxes); % Check if it is a Axes handle. For old Matlab versions, the test with isnan is used.
    haxes = gca;
    hfig = get(haxes,'Parent');
else
    haxes = data.options.haxes;
    hfig = get(haxes,'Parent');
end

data.options.hfig = hfig;
data.options.haxes = haxes; % Save axes handle for post-processing of the plot

% normalize fs
fixscale = data.flags.fixscale;
if normfs
    fscale = 1 / fs;
elseif fs>10e3&&fs<10e6&&~fixscale
    fscale = 1/1e3;
elseif fs>10e6&&~fixscale
    fscale = 1/1e6;
else
    fscale = 1;
end

% Hold on
set(haxes,'Next','Add');


if strcmpi(scale,'sindb')
    % Spectrum without signal
    plot(haxes,f*fscale,10*log10(psd*dfbin/PowerFS),'LineStyle','-','color',c_psd);
    % Plot accumulated IBN
    if data.options.HighlightIB
        plot(haxes,f(bins_fband)*fscale,10*log10(psd_IB_accumulated),'LineStyle','-','color',c_ibn);
    end
    % Highlight signal
    if data.options.HighlightSignal
        for ni = 1:Nsig
            plot(haxes,f(bins_signal{ni})*fscale,10*log10(psd(bins_signal{ni})*dfbin/PowerFS),'Marker','x','LineStyle','-','color',c_signal);
        end
    end
    grid on;
    ylabel('PSD (dBFS/bin)');
end

if strcmpi(scale,'sin')
    % Spectrum without signal *NG*fbin/CG^2
    scalePower = NG*dfbin/CG^2*PowerFS*PowerRescaleFS;
    plot(haxes,f*fscale,sqrt(2*psd*scalePower),'LineStyle','-','color',c_psd);
    % Plot accumulated IBN
    if data.options.HighlightIB
        plot(haxes,f(bins_fband)*fscale,sqrt(2*psd_IB_accumulated),'LineStyle','-','color',c_ibn);
    end
    % Highlight signal
    if data.options.HighlightSignal
        for ni = 1:Nsig
            plot(haxes,f(bins_signal{ni})*fscale,sqrt(2*psd(bins_signal{ni})*scalePower),'Marker','x','LineStyle','-','color',c_signal);
        end
    end
    grid(haxes,'on');
    ylabel(haxes,['ASD (' unit_str ')']);
    set(haxes,'YScale','log');
end

if strcmpi(scale,'nbwdb')
    % Spectrum without signal
    plot(haxes,f*fscale,10*log10(psd*NG/PowerFS),'LineStyle','-','color',c_psd);
    % Plot accumulated IBN
    if data.options.HighlightIB
        plot(haxes,f(bins_fband)*fscale,10*log10(psd_IB_accumulated),'LineStyle','-','color',c_ibn);
    end
    % Highlight signal
    if data.options.HighlightSignal
        for ni = 1:Nsig
            plot(haxes,f(bins_signal{ni})*fscale,10*log10(psd(bins_signal{ni})*NG/PowerFS),'Marker','x','LineStyle','-','color',c_signal);
        end
    end
    grid on;
    ylabel('PSD (dBFS/NBW)');
end

if strcmpi(scale,'nbw')
    Wnorm = 1;% length(W)*(W'*W)/(norm(W,1)^2);
    % Wnorm = length(W)/(norm(W,1)^2);
    % Spectrum normalized (dBFS/NBW)  without signal
    if show_IQ
        psd_I = data.spectrum.psd_I;
        psd_Q = data.spectrum.psd_Q;
        f_IQ = data.spectrum.f_IQ;
        plot(haxes,f_IQ*fscale,(psd_I/PowerFS*PowerRescaleFS*Wnorm).^0.5,'LineStyle','-','color',c_psd_I,'DisplayName','Spectrum_I');
        hold on;
        plot(haxes,f_IQ*fscale,(psd_Q/PowerFS*PowerRescaleFS*Wnorm).^0.5,'LineStyle','-','color',c_psd_Q,'DisplayName','Spectrum_Q');
    end
    plot(haxes,f*fscale,(psd/PowerFS*PowerRescaleFS*Wnorm).^0.5,'LineStyle','-','color',c_psd,'DisplayName','Spectrum');
    % Accumulated noise
    if data.options.HighlightIB
        plot(haxes,f(bins_fband)*fscale,(psd_IB_accumulated/PowerFS*PowerRescaleFS*Wnorm).^0.5,c_ibn,'DisplayName','Integrated IBN');
    end
    % Highlight signal peak
    if data.options.HighlightSignal
        for ni = 1:Nsig
            loglog(haxes,f(bins_signal{ni})*fscale,(psd(bins_signal{ni})/PowerFS*PowerRescaleFS*Wnorm).^0.5,'Marker','x','LineStyle','-','color',c_signal,'DisplayName','Signal');
        end
    end
    grid(haxes,'on');
    ylabel(haxes,['PSD (' unit_str '/rtHz)']);
    set(haxes,'YScale','log');
end

if strcmpi(scale,'spsd')
    color_spcc = data.color.parameters.color_spcc;
    set(haxes,'NextPlot','Add');
    plot(haxes,f*fscale,real(psd),'linestyle','none','color',color_spcc,'Marker','x','DisplayName','Rate space');
    plot(haxes,f*fscale,imag(psd),'linestyle','none','color',color_spcc,'Marker','o','DisplayName','Quad space');
    grid(haxes,'on');
    axis(haxes,'tight');
    ylabel('Spectral Cross-Correlation (-)');
end


% Plot labels & text
if normfs
    xlabel(haxes,'Frequency (f_s)');
    xlim(haxes,[1/Np 1/2]*fscale);
    fs = 1; %#ok<NASGU>
elseif fs>10e3&&fs<10e6&&~fixscale
    xlabel(haxes,'Frequency (kHz)');
    xlim(haxes,[dfbin fs/2]*fscale);
elseif fs>10e6&&~fixscale
    xlabel(haxes,'Frequency (MHz)');
    xlim(haxes,[dfbin fs/2]*fscale);
else
    xlabel(haxes,'Frequency (Hz)');
    xlim(haxes,[dfbin fs/2]*fscale);
end

title(haxes,str_title,'Interpreter','none')

stats_position = data.options.stats_position;
switch lower(stats_position)
    case 'north'
        xpos = 0.5; ypos = 0.9; hor_str = 'center'; ver_str = 'top'; % hor: {left} | center | right; ver: top | cap | {middle} | baseline | bottom
    case 'south'
        xpos = 0.5; ypos = 0.1; hor_str = 'center'; ver_str = 'bottom';
    case 'west'
        xpos = 0.1; ypos = 0.5; hor_str = 'left'; ver_str = 'middle';
    case 'east'
        xpos = 0.9; ypos = 0.1; hor_str = 'right'; ver_str = 'middle';
    case 'northwest'
        xpos = 0.1; ypos = 0.9; hor_str = 'left'; ver_str = 'top';
    case 'southwest'
        xpos = 0.1; ypos = 0.1; hor_str = 'left'; ver_str = 'bottom';
    case 'southeast'
        xpos = 0.9; ypos = 0.1; hor_str = 'right'; ver_str = 'bottom';
    case 'northeast'
        xpos = 0.9; ypos = 0.9; hor_str = 'right'; ver_str = 'top';
    otherwise
        error('plotFunction:PositionStatField','Position %s not supported',stats_position)
end

if stats
    if strcmpi(scale,'nbw')||strcmpi(scale,'sin')
        sig_amp = sqrt(data.signal.P_signal*2/PowerFS*PowerRescaleFS);
        sel = fsig==0;
        sig_amp(sel) = sqrt(data.signal.P_signal(sel)/PowerFS*PowerRescaleFS);
        data.signal.Amplitude = sig_amp;
        noise_rms = sqrt(data.ibn.P_IBN/PowerFS*PowerRescaleFS);
        pvalue = data.spectrum.parameters.pvalue;
        
        % Noise density within e.g. 1...20Hz range
        fmin = data.spectrum.AverageNoiseDensity.fmin;
        fmax = data.spectrum.AverageNoiseDensity.fmax;
        Smean = data.spectrum.AverageNoiseDensity.AND;
        [scale_noise_rms, str_noise_rms] = FncScaleSignal(noise_rms);
        [scale_sig_amp, str_sig_amp] = FncScaleSignal(sig_amp(1));
        if data.flags.calculate_IQ
            sig_amp_I =  sqrt(data.signal.P_signal_I_renorm); % Signal is at DC. Omit factor of 2
            sig_amp_Q =  sqrt(data.signal.P_signal_Q_renorm); % Signal is at DC. Omit factor of 2
            noise_rms_I = sqrt(data.ibn.P_IBN_I_renorm);
            noise_rms_Q = sqrt(data.ibn.P_IBN_Q_renorm);
            phi_demod = data.timesignal.parameters.phi_demod;
%             [scale_noise_rms_I, str_noise_rms_I] = FncScaleSignal(noise_rms_I);
%             [scale_noise_rms_Q, str_noise_rms_Q] = FncScaleSignal(noise_rms_Q);
            scale_noise_rms_I = scale_noise_rms;
            scale_noise_rms_Q = scale_noise_rms;
            str_noise_rms_I = str_noise_rms;
            str_noise_rms_Q = str_noise_rms;
%             [scale_sig_amp_I, str_sig_amp_I] = FncScaleSignal(sig_amp_I);
%             [scale_sig_amp_Q, str_sig_amp_Q] = FncScaleSignal(sig_amp_Q);
            scale_sig_amp_I = scale_sig_amp;
            scale_sig_amp_Q = scale_sig_amp;
            str_sig_amp_I = str_sig_amp;
            str_sig_amp_Q = str_sig_amp;
            
            if pvalue>0
                Neff = size(psd_IB_accumulated,1);
                ConfInt = [1 sqrt([(Neff-1)/chi2inv(1-pvalue/2,Neff-1) (Neff-1)/chi2inv(pvalue/2,Neff-1)])];
                data.ibn.ConfInt = ConfInt;
                text('Parent',haxes,'units','normalized','HorizontalAlignment',hor_str,...
                    'VerticalAlignment',ver_str,'Position',[xpos ypos],'String',...
                    {sprintf('Confidence level of 1-%1.2f%%:',pvalue*100),...
                    sprintf('IBN = %3.4f [%3.4f %3.4f] (%s,rms)',ConfInt*noise_rms*scale_noise_rms,[str_noise_rms unit_str]),...
                    sprintf('IBN_I = %3.4f [%3.4f %3.4f] (%s,rms)',ConfInt*noise_rms_I*scale_noise_rms_I,[str_noise_rms_I unit_str]),...
                    sprintf('IBN_Q = %3.4f [%3.4f %3.4f] (%s,rms)',ConfInt*noise_rms_Q*scale_noise_rms_Q,[str_noise_rms_Q unit_str]),...
                    sprintf('A_{sig} = %3.3f (%s)',sig_amp(1)*scale_sig_amp,[str_sig_amp unit_str]),...
                    sprintf('A_{sig,I} = %3.3f (%s)',sig_amp_I*scale_sig_amp_I,[str_sig_amp_I unit_str]),...
                    sprintf('A_{sig,Q} = %3.3f (%s)',sig_amp_Q*scale_sig_amp_Q,[str_sig_amp_Q unit_str]),...
                    sprintf('Phi_{demod} = %3.3f (°)',phi_demod*180/pi),...
                    sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean,unit_str)},...
                    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5,'Parent',haxes,'Tag','StatisticText');
            else
                text('Parent',haxes,'units','normalized','HorizontalAlignment',hor_str,...
                    'VerticalAlignment',ver_str,'Position',[xpos ypos],'String',...
                    {sprintf('IBN = %3.4f (%s,rms)',noise_rms*scale_noise_rms,[str_noise_rms unit_str]),...
                    sprintf('IBN_I = %3.4f (%s,rms)',noise_rms_I*scale_noise_rms_I,[str_noise_rms_I unit_str]),...
                    sprintf('IBN_Q = %3.4f (%s,rms)',noise_rms_Q*scale_noise_rms_Q,[str_noise_rms_Q unit_str]),...
                    sprintf('A_{sig} = %3.3f (%s)',sig_amp(1)*scale_sig_amp,[str_sig_amp unit_str]),...
                    sprintf('A_{sig,I} = %3.3f (%s)',sig_amp_I*scale_sig_amp_I,[str_sig_amp_I unit_str]),...
                    sprintf('A_{sig,Q} = %3.3f (%s)',sig_amp_Q*scale_sig_amp_Q,[str_sig_amp_Q unit_str]),...
                    sprintf('Phi_{demod} = %3.3f (°)',phi_demod*180/pi),...
                    sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean,unit_str)},...
                    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5,'Parent',haxes,'Tag','StatisticText');
            end
        else
            if pvalue>0
                Neff = size(psd_IB_accumulated,1);
                ConfInt = [1 sqrt([(Neff-1)/chi2inv(1-pvalue/2,Neff-1) (Neff-1)/chi2inv(pvalue/2,Neff-1)])];
                data.ibn.ConfInt = ConfInt;
                text('Parent',haxes,'units','normalized','HorizontalAlignment',hor_str,...
                    'VerticalAlignment',ver_str,'Position',[xpos ypos],'String',...
                    {sprintf('Confidence level of 1-%1.2f%%:',pvalue*100),...
                    sprintf('IBN = %3.4f [%3.4f %3.4f] (%s,rms)',ConfInt*noise_rms*scale_noise_rms,[str_noise_rms unit_str]),...
                    sprintf('A_{sig} = %3.3f (%s)',sig_amp(1)*scale_sig_amp,[str_sig_amp unit_str]),...
                    sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean*1000,['m' unit_str])},...
                    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5,'Parent',haxes,'Tag','StatisticText');
%                     sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean*scale_noise_rms,[str_noise_rms unit_str])},...
            else
                text('Parent',haxes,'units','normalized','HorizontalAlignment',hor_str,...
                    'VerticalAlignment',ver_str,'Position',[xpos ypos],'String',...
                    {sprintf('IBN = %3.4f (%s,rms)',noise_rms*scale_noise_rms,[str_noise_rms unit_str]),...
                    sprintf('A_{sig} = %3.3f (%s)',sig_amp(1)*scale_sig_amp,[str_sig_amp unit_str]),...
                    sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean*1000,['m' unit_str])},...
                    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5,'Parent',haxes,'Tag','StatisticText');
                %                     sprintf('Noise [%1.0f-%1.0fHz] = %3.3f (%s/rtHz)',fmin,fmax,Smean*scale_noise_rms,[str_noise_rms unit_str])},...
            end
        end
    else
        text('Parent',haxes,'units','normalized','Position',[xpos ypos],'HorizontalAlignment',hor_str,...
            'VerticalAlignment',ver_str,'String',{['IBN = ',num2str(data.ibn.P_IBN_dBFS,'%6.2f'),' (dBFS);'],...
            ['P_{sig} = ',num2str(data.signal.P_signal_dBFS,'%6.2f'),' (dBFS)']},...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'Parent',haxes,'Tag','StatisticText');
    end
end

if flabel
    if (segment_no == 1)
        text('Parent',haxes,'units','normalized','Position',[xpos ypos],'HorizontalAlignment',hor_str,...
            'VerticalAlignment',ver_str,'String', {'Blackman-Harris window',...
            ['2\^' num2str(floor(log2(Np))) ' pt FFT'], [num2str(segment_no) ' segment, no overlap' ]},...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5,'Parent',haxes);
    else
        text('Parent',haxes,'units','normalized','Position',[xpos ypos],'HorizontalAlignment',hor_str,...
            'VerticalAlignment',ver_str,'String', {'Blackman-Harris window',...
            ['2\^' num2str(floor(log2(Np))) ' pt FFT'],...
            [num2str(segment_no) ' segments, ' num2str(overlap_string) '% overlap']},...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'Parent',haxes);
    end
end

function data = ChkArgs(args,ni)
% Default structure
data.timesignal = struct;
data.timesignal.parameters = struct;
data.spectrum = struct;
data.spectrum.parameters = struct;
data.options = struct;
data.flags = struct;
data.info = struct;
data.color = struct;
data.color.parameters = struct;

data.flags.forcedirectplot = false;
data.flags.init = false;

if ni == 0
    error('Syntax error:at least one argument required. Please type ''help plotFunction'' into the command window.');
elseif isa(args{1},'struct')
    data = args{1};
    data.flags.forcedirectplot = true;
    return;
elseif isa(args{1},'char')&&strcmpi(args{1},'init')
    data.flags.init = true;
elseif ~isa(args{1},'numeric')
    error('Syntax error:First argument is not numeric')
elseif ~isvector(args{1})
    error('Syntax error::First argument is not an array')
else
    v = args{1};
    % Check if vector is column or row
    [Np,n2] = size(v);
    if Np~=1&&n2~=1
        error('Syntax error:signal is multidimensional')
    end
    
    if Np<n2
        v = v';
    end
    data.timesignal.v = v;
end

if ni==2
    if isa(args{2},'struct')
        data = args{2};
        data.timesignal.v = args{1};
    else
        error('Syntax error:command called with two arguments, where non of them is a structure')
    end
else
    % Check syntax {yout,'str_1',var_1,'str_2',var_2,...}
    for i = 2:2:ni-1
        if ~isa(args{i},'char')
            error('Syntax error::string assumed')
        end
    end
    
    data = ChkString(data,args,'fs','timesignal','numeric','sample frequency is not numeric');
    data = ChkString(data,args,'SignalSizeCorrection','timesignal','char','SignalSizeCorrection is not a string');
    data = ChkString(data,args,'SignalSizeCorrectionMethod','timesignal','char','SignalSizeCorrectionMethod is not a string');
    data = ChkString(data,args,'normfs','flags','logical','norm fs is not logical');
    data = ChkString(data,args,'scale','options','char','PSD scaling is not a string');
    data = ChkString(data,args,'HighlightSignal','options','logical','HighlightSignal is noot logical');
    data = ChkString(data,args,'HighlightIB','options','logical','HighlightIB is noot logical');
    data = ChkString(data,args,'finephaseno','options','numeric','finephaseno is not a number.');
    data = ChkString(data,args,'stats_position','options','char','Position of the statistic text is not a string.');
    data = ChkString(data,args,'Np','timesignal','numeric','number of points is not an integer');
    data = ChkString(data,args,'segment_no','spectrum','numeric','segment number is not an integer');
    data = ChkString(data,args,'pvalue','spectrum','numeric','pvalue is not a number');
    data = ChkString(data,args,'iq_noise','spectrum','char','iq_noise is not a character');
    data = ChkString(data,args,'OSR','timesignal','numeric','OSR is not numeric');
    data = ChkString(data,args,'fb','timesignal','numeric','fb is not numeric');
    data = ChkString(data,args,'fsig','timesignal','numeric','signal frequency is not a number');
    data = ChkString(data,args,'FS','spectrum','numeric','full scale is not numeric');
    data = ChkString(data,args,'killbins','spectrum','numeric','killbins is not ');
    data = ChkString(data,args,'RescaleFS','spectrum','numeric','Rescaled full scale is not numeric');
    data = ChkString(data,args,'new_fig','flags','logical','boolean expression assumed');
    data = ChkString(data,args,'plot_fft','flags','logical','boolean expression assumed');
    data = ChkString(data,args,'haxes','options','numeric','Handle for axes assumed','matlab.graphics.axis.Axes');
    data = ChkString(data,args,'f0','timesignal','numeric','center frequency is not numeric');
    data = ChkString(data,args,'stats','flags','logical','stats is not logical');
    data = ChkString(data,args,'flabel','flags','logical','flabel is not logical');
    data = ChkString(data,args,'ReducedMemoryStorage','flags','logical','ReducedMemoryStorage is not logical');
    data = ChkString(data,args,'DigitalCorrection','flags','logical','DigitalCorrection is not logical');
    data = ChkString(data,args,'fmod','spectrum','numeric','fmod is not numeric');
    data = ChkString(data,args,'no_nz_bins','spectrum','numeric','bins is not numeric');
    data = ChkString(data,args,'verbose','flags','logical','verbose is not true or false');
    data = ChkString(data,args,'show_IQ','flags','logical','show_IQ is not true or false');
    data = ChkString(data,args,'fixscale','flags','logical','fixscale is not true or false');
    data = ChkString(data,args,'calculate_IQ','flags','logical','calculate_IQ is not true or false');
    data = ChkString(data,args,'phi_demod','timesignal','numeric','phi_demod is not numeric','char');
    data = ChkString(data,args,'predict_phi_demod','options','char','predict_phi_demod is not I or Q');
    data = ChkString(data,args,'SpectralCrossCorrelation','flags','logical','SpectralCrossCorrelation is not true of false');
    data = ChkString(data,args,'FocusInBand','options','numeric','FocusInBand is not true of false');
    data = ChkString(data,args,'color_spcc','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'color_signal','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'color_ibn','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'color_spectrum','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'color_spectrum_rate','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'color_spectrum_quad','color','numeric','Color vector is not numeric');
    data = ChkString(data,args,'unit','timesignal','char','Unit is not a string');
    data = ChkString(data,args,'CorrelationFunction','options','char','String for CorrelationFunction assumed');
    data = ChkString(data,args,'str_title','options','char','Title must be a string');
    data = ChkString(data,args,'MaximumFrequency','options','numeric','MaximumFreuqnecy is not numeric.');
    data = ChkString(data,args,'FrequencyResolution','options','numeric','FrequencyResolution is not numeric.');
end

data = ChkStruct(data);


function data = ChkString(data,args,string,str_group,type,err_msg,type_2)
if nargin < 7
    type_2 = type;
end

post = find(strcmp(args(2:2:end),string),1);
pos = 2*post;
if isempty(pos)
elseif ~(isa(args{pos+1},type) || isa(args{pos+1},type_2))
    error(['Syntax error::',err_msg])
elseif (isa(args{pos+1},type) && strcmp('numeric',type) && (any(isnan(args{pos+1})) || any(isinf(args{pos+1}))))
    error(['Syntax error::',err_msg])
else
    if strcmpi(str_group,'options')||strcmpi(str_group,'flags')
        data.(str_group).(string)=args{pos+1};
    else
        data.(str_group).parameters.(string)=args{pos+1};
    end
end


function data = ChkStruct(data)
data.info.txt = '\n';
if ~isfield(data.flags,'normfs')
    data.flags.normfs = 0;
end

if ~isfield(data.options,'HighlightSignal')
    data.options.HighlightSignal = true;
end

if ~isfield(data.options,'HighlightIB')
    data.options.HighlightIB = true;
end

if ~isfield(data.options,'scale')
    data.options.scale = 'nbw';
end

if ~isfield(data.options,'predict_phi_demod')
    data.options.predict_phi_demod = 'I';
    data.flags.flag_predict_phi_demod = false;
else
    data.flags.flag_predict_phi_demod = true;
end

if ~isfield(data.timesignal.parameters,'phi_demod')
    data.timesignal.parameters.phi_demod = 0;
elseif ischar(data.timesignal.parameters.phi_demod)
    data.options.predict_phi_demod = data.timesignal.parameters.phi_demod;
    data.flags.flag_predict_phi_demod = true;
    data.timesignal.parameters.phi_demod = 0;
end

if ~isfield(data.flags,'show_IQ')
    data.flags.show_IQ = false;
elseif data.flags.show_IQ
    data.flags.calculate_IQ = true;
end

if ~isfield(data.flags,'calculate_IQ')
    data.flags.calculate_IQ = false;
end

if ~isfield(data.flags,'ReducedMemoryStorage')
    data.flags.ReducedMemoryStorage = false;
end

if ~isfield(data.timesignal.parameters,'Np')
    data.timesignal.parameters.Np = length(data.timesignal.v);
end

if ~isfield(data.timesignal.parameters,'SignalSizeCorrection')
    data.timesignal.parameters.SignalSizeCorrection = 'ForceTruncationPower2';
end

if ~isfield(data.timesignal.parameters,'SignalSizeCorrectionMethod')
    data.timesignal.parameters.SignalSizeCorrectionMethod = 'zeros';
end

if ~isfield(data.spectrum.parameters,'segment_no')
    data.spectrum.parameters.segment_no = 1;
elseif mod(data.spectrum.parameters.segment_no,2)
    data.info.txt=[data.info.txt,'Warning:signal is not located on a bin\n'];
end

if ~isfield(data.timesignal.parameters,'fs')
    data.timesignal.parameters.fs = 1;
end

if ~isfield(data.timesignal.parameters,'fb')
    data.timesignal.parameters.fb = data.timesignal.parameters.fs/4000;
end

if ~isfield(data.timesignal.parameters,'fsig')
    data.timesignal.parameters.fsig = data.timesignal.parameters.fs/(16);
end

if ~isfield(data.spectrum.parameters,'FS')
    data.spectrum.parameters.FS = 2;
end

if ~isfield(data.timesignal.parameters,'unit')
    % data.timesignal.parameters.unit = 'dps';
    data.timesignal.parameters.unit = 'AU';
end

if ~isfield(data.spectrum.parameters,'RescaleFS')
    data.spectrum.parameters.RescaleFS = 2;
end

if ~isfield(data.spectrum.parameters,'killbins')
    data.spectrum.parameters.killbins = [];
end

if ~isfield(data.spectrum.parameters,'overlapRatio')
    data.spectrum.parameters.overlapRatio = 0.5;
end

if ~isfield(data.flags,'plot_fft')
    data.flags.plot_fft = true;
end

if ~isfield(data.spectrum.parameters,'pvalue')
    data.spectrum.parameters.pvalue = 0.01;
end

if ~isfield(data.flags,'new_fig')
    data.flags.new_fig = false;
elseif data.flags.new_fig
    data.flags.plot_fft = true;
end

if ~isfield(data.options,'haxes')
    data.options.haxes = NaN;
end

if ~isfield(data.spectrum.parameters,'no_nz_bins')
    data.spectrum.parameters.no_nz_bins = 3;
end

if ~isfield(data.spectrum.parameters,'iq_noise')
    data.spectrum.parameters.iq_noise = 'two_sided';
end

if ~isfield(data.flags,'flabel')
    data.flags.flabel = false;
end

if ~isfield(data.options,'str_title')
    data.options.str_title = '';
end

if ~isfield(data.options,'CorrelationFunction')
    data.options.CorrelationFunction = 'cxcorr';
end

if ~isfield(data.options,'MaximumFrequency')
    data.options.MaximumFrequency = 5*25e3;
end

if ~isfield(data.options,'FrequencyResoulution')
    data.options.FrequencyResoulution = 1;
end

if ~isfield(data.timesignal.parameters,'f0')
    data.timesignal.parameters.f0 = 0;
end

if isfield(data.timesignal.parameters,'OSR')
    fprintf('OSR will be obsolete in a future version of plotFunction.\n');
    if data.timesignal.parameters.f0 == 0
        data.timesignal.parameters.fb = data.timesignal.parameters.fs/(2*data.timesignal.parameters.OSR);
    else
        data.timesignal.parameters.fb = data.timesignal.parameters.fs/(4*data.timesignal.parameters.OSR);
    end
end

if ~isfield(data.flags,'stats')
    data.flags.stats = false;
end

if ~isfield(data.options,'stats_position')
    data.options.stats_position = 'southwest';
end

if ~isfield(data.options,'finephaseno')
    data.options.finephaseno = 10;
end

if ~isfield(data.flags,'fixscale')
    data.flags.fixscale = false;
end

if ~isfield(data.options,'verbose')
    data.info.verbose = false;
else
    data.info.verbose = data.options.verbose;
    data = rmfield(data,'verbose');
end

if ~isfield(data.flags,'SpectralCrossCorrelation')
    data.flags.SpectralCrossCorrelation = false;
elseif data.flags.SpectralCrossCorrelation
    data.flags.stats=false;
    data.flags.plot_fft = true;
    data.options.scale = 'spsd';
end

if ~isfield(data.options,'FocusInBand')
    data.options.FocusInBand = 0;
else
    if data.options.FocusInBand < 0
        data.options.FocusInBand = abs(data.options.FocusInBand);
    end
end

if ~isfield(data.color.parameters,'color_spcc')
    data.color.parameters.color_spcc = [1 0 0];
end

if ~isfield(data.spectrum.parameters,'Nsig')
    data.spectrum.parameters.Nsig = 1;
end

if ~isfield(data.flags,'DigitalCorrection')
    data.flags.DigitalCorrection = false;
end

if ~isfield(data.spectrum.parameters,'fmod')
    f0 = data.timesignal.parameters.f0;
    data.spectrum.parameters.fmod = [f0 2*f0 3*f0 4*f0];
end
%data.color.psd = 'k'; %'r';
if ~isfield(data.color.parameters,'color_spectrum')
    data.color.psd = 'k';
else
    data.color.psd = data.color.parameters.color_spectrum;
end


%data.color.psd_I = 'r'; %'r';
if ~isfield(data.color.parameters,'color_spectrum_rate')
    data.color.psd_I = 'r';
else
    data.color.psd_I = data.color.parameters.color_spectrum_rate;
end

%data.color.psd_Q = 'm'; %'r';
if ~isfield(data.color.parameters,'color_spectrum_quad')
    data.color.psd_Q = 'm';
else
    data.color.psd_Q = data.color.parameters.color_spectrum_quad;
end

%data.color.signal = 'g'; %'g*-';
if ~isfield(data.color.parameters,'color_signal')
    data.color.signal = 'g';
else
    data.color.signal = data.color.parameters.color_signal;
end

data.color.ibn = 'b';
data.spectrum.PowerFS = (data.spectrum.parameters.FS/2)^2; % Power of a sysmetric rect signal. Sine wave has half of that!
data.spectrum.PowerRescaleFS = (data.spectrum.parameters.RescaleFS/2)^2;

function data = FncCalculateBINs(data)
fsigs = abs(data.timesignal.parameters.fsig-data.timesignal.parameters.fs*round(data.timesignal.parameters.fsig/data.timesignal.parameters.fs));
bin_signal = floor(fsigs/data.spectrum.dfbin)+1;
bin_f0 = floor(data.timesignal.parameters.f0/data.spectrum.dfbin)+1;
bin_fband = floor(data.timesignal.parameters.fb/data.spectrum.dfbin)+1;

newfb = data.spectrum.dfbin * bin_fband;
fb = data.timesignal.parameters.fb;
if abs((newfb/fb-1))>0.05 || abs(newfb-fb)>5
    warning('plotFunction:ChkBandwidth','Due to frequency resolution of the spectrum, the effective bandwidth is %1.2f',newfb)
end

if (bin_fband)<10
    bin_fband = 10;
    data.timesignal.parameters.fb = data.spectrum.dfbin*bin_fband;
    warning('plotFunction:ChkNpIBN','band-width is increased to %2.1f Hz.',data.timesignal.parameters.fb);
end

Nsig = length(bin_signal);
bin_signal_min = max((bin_signal-data.spectrum.parameters.no_nz_bins),1);
bin_signal_max = min((bin_signal+data.spectrum.parameters.no_nz_bins),floor(data.timesignal.parameters.Np/2));
bins_signal_t = cell(Nsig,1);
parfor i=1:Nsig
    bins_signal_t{i} = bin_signal_min(i):bin_signal_max(i);
end
bins_signal = unique([bins_signal_t{:}]);
bins_fband_t = max((bin_f0-bin_fband),1):min((bin_f0+bin_fband),floor(data.timesignal.parameters.Np/2));
killbins = data.spectrum.parameters.killbins;
killbins(killbins>floor(data.timesignal.parameters.Np/2))=[];
bins_fband_t(killbins) = [];
bins_fband = setdiff(bins_fband_t,bins_signal);
bins_fband_intersect = intersect(bins_fband_t,bins_signal);

data.spectrum.parameters.bins_signal = bins_signal_t;
data.spectrum.parameters.Nsig = Nsig;
data.spectrum.parameters.bins_fband = bins_fband;
data.spectrum.parameters.bins_fband_intersect = bins_fband_intersect;
