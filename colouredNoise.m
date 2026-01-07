function [t,y] = colouredNoise(tL,fS,fASD,ASD,seed,fsin,Asin,phisin,plot_flag)
t   = [0:1/fS:tL];           % time vector in s
N   = length(t);             % number of sampling points
N1  = floor(N/2)+1;
f   = [0:fS/N:fS*(1-1/N)];   % frequency vector in Hz
flg = log10(f);

%% check input parameters
if fASD(1) == 0
   fASD(1) = f(2);
end
if min(fASD) < 0
   error('fASD < 0')
end
if (min(ASD) <= 0) && (max(ASD) > 0)
   error('ASD <= 0')
end   

%% calculation of noise
if sum(ASD) == 0
    y = zeros(1,N);
else
    fASDlg  = log10(fASD);
    DfASDlg = diff(fASDlg);
    ASDlg   = log10(ASD);
    DASDlg  = diff(ASDlg);
    Ades    = ones(1,N1).*(min(ASDlg)-16);
    nlast   = min([find(f >= fASD(1),1,'first') N1]);
    for k = 1:length(DASDlg)
       n             = min([find(f >= fASD(k+1),1,'first') N1]);
       Ades(nlast:n) = ASDlg(k)+DASDlg(k)/DfASDlg(k)*(flg(nlast:n)-fASDlg(k));
       nlast         = n;
    end
    Ades    = 10.^Ades;
    Ades(1) = 0; % Ades(2);
    rng(seed)
    phi     = [0 (rand(1,N1-2)-.5).*(2*pi) 0];
    A       = Ades.*exp(1i*phi);
    A       = [A fliplr(conj(A(2:ceil(N/2))))]./sqrt(2);
    y       = real(N.*ifft(A.*sqrt(fS/N)));
end

%% add additional tones
y = y+(sum((ones(N,1)*Asin).*sin(mod(t'*((2*pi).*fsin)+ones(N,1)*(phisin.*(pi/180)),2*pi)),2))';

%% plot
if plot_flag == 1
   Y   = fft(y)./N;
   PSD = (N/fS).*abs(Y).^2;
   figure; plot(t,y); grid on; xlabel('t / s'); ylabel('y');
   figure; loglog(f,sqrt(2.*PSD)); grid on; xlim([0 fS/2]);
   xlabel('f / Hz'); ylabel('single sideband ASD / 1/rtHz');
end
