function varargout = spectrokqgain(V,Y,varargin)
%
% spectrokqgain(X,Y,NOPOINTS,NOVERLAP,Ts)
% [kq] = spectrohist(X,Y,NOPOINTS,NOVERLAP,Ts)
% Inputs:
% V ------------------------------------------- output signal of quantizer
% Y ------------------------------------------- Input signal of quantizer
% NOPOINTS ------------------------------------ Separation into segments that have NOPOINTS elements [2^9]
% NOVERLAP ------------------------------------ Number of elements that overlap between consecutive segements [2^8]
% Ts ------------------------------------------ Sampling time [1]

%
% Outputs:
% kq ------------------------------------------ Quantizer gain over time

% AE/ESE3-Buhmann

switch nargin,
    case 0,
        error('spectrokqgain:chk_args','At least two input singals are required.')
    case 1,
        error('spectrokqgain:chk_args','At least two input singals are required.')
    case 2,
        NOPOINTS = 2^9;
        NOVERLAP = 2^8;
        Ts = 1;
    case 3,
        NOPOINTS = varargin{1};
        NOVERLAP = 2^8;
        Ts = 1;
    case 4,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        Ts = 1;
    case 5,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        Ts = varargin{3};
end

nx = length(V);
ncol = fix((nx-NOVERLAP)/(NOPOINTS-NOVERLAP));
colindex = 1 + (0:(ncol-1))*(NOPOINTS-NOVERLAP);
rowindex = (1:NOPOINTS)';
Vt = zeros(NOPOINTS,ncol);
Vt(:) = V(rowindex(:,ones(1,ncol))+colindex(ones(NOPOINTS,1),:)-1);
Yt = zeros(NOPOINTS,ncol);
Yt(:) = Y(rowindex(:,ones(1,ncol))+colindex(ones(NOPOINTS,1),:)-1);

kq = zeros(1,ncol);

parfor ni = 1:ncol
    kq(1,ni) = (Yt(:,ni)'*Vt(:,ni))/((Yt(:,ni)'*Yt(:,ni)));
end

switch nargout,
    case 0,
        plot(Ts*(NOPOINTS/2)*(1:ncol),kq,'k-');
        axis xy; axis tight;
        ylabel('Quantizer gain (-)');
        xlabel('Time')
    case 1,
        varargout = {kq};
    case 3,
        varargout = {kq,NOPOINTS,NOVERLAP};
    case 4,
        varargout = {kq,NOPOINTS,NOVERLAP,ncol};
    case 5,
        varargout = {kq,NOPOINTS,NOVERLAP,ncol,Ts};
end
