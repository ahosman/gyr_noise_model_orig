function varargout = spectrohist(X,varargin)
%
% spectrohist(X,NOPOINTS,NOVERLAP,NOBIN,Ts,irange)
% [n,bins] = spectrohist(X,NOPOINTS,NOVERLAP,NOBIN,Ts,irange)
% Inputs:
% X ------------------------------------------- Input signal
% NOPOINTS ------------------------------------ Separation into segments that have NOPOINTS elements [2^9]
% NOVERLAP ------------------------------------ Number of elements that overlap between consecutive segements [2^8]
% NOBIN --------------------------------------- Number of bins of the histogram [100]
% Ts ------------------------------------------ Sampling time [1]
% irange -------------------------------------- Range of the input signal [-1 1]

%
% Outputs:
% n ------------------------------------------- number of elements inside bin
% bins ---------------------------------------- bin center

% AE/ESE3-Buhmann

switch nargin,
    case 0,
        error('spectrohist:chk_args','At least a input singal is required.')
    case 1,
        NOPOINTS = 2^9;
        NOVERLAP = 2^8;
        NOBIN = 100;
        Ts = 1;
        irange = [-1 1];
    case 2,
        NOPOINTS = varargin{1};
        NOVERLAP = 2^8;
        NOBIN = 100;
        Ts = 1;
        irange = [-1 1];
    case 3,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        NOBIN = 100;
        Ts = 1;
        irange = [-1 1];
    case 4,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        NOBIN = varargin{3};
        Ts = 1;
        irange = [-1 1];
    case 5,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        NOBIN = varargin{3};
        Ts = varargin{4};
        irange = [-1 1];
    case 6,
        NOPOINTS = varargin{1};
        NOVERLAP = varargin{2};
        NOBIN = varargin{3};
        Ts = varargin{4};
        irange = varargin{5};
end

nx = length(X);
ncol = fix((nx-NOVERLAP)/(NOPOINTS-NOVERLAP));
colindex = 1 + (0:(ncol-1))*(NOPOINTS-NOVERLAP);
rowindex = (1:NOPOINTS)';
xin = zeros(NOPOINTS,ncol);
xin(:) = X(rowindex(:,ones(1,ncol))+colindex(ones(NOPOINTS,1),:)-1);

NOE = zeros(NOBIN,ncol);
BINV = zeros(NOBIN,ncol);

BINS = linspace(irange(1),irange(2),NOBIN);

for ni = 1:ncol
    [NOE(:,ni) BINV(:,ni)] = hist(xin(:,ni),BINS);
%     NOE(NOE(:,ni)==0,ni)=NaN;
    norm_pdf = [diff(BINV(:,ni)') 0]*NOE(:,ni);
    NOE(:,ni) = NOE(:,ni)/norm_pdf;
end

switch nargout,
    case 0,
        [X,Y] = meshgrid(1:ncol,1:NOBIN);
        surf(X*Ts,BINV,NOE,'EdgeColor','none');
        axis xy; axis tight;
        colormap(jet);
        view(0,90);
        ylabel('Input signal');
        xlabel('Time')
        zlabel('Probability')
%         colorbar('location','eastoutside')
        colorbar('location','east')
    case 1,
        varargout = {NOE};
    case 2,
        varargout = {NOE,BINV};
    case 3,
        varargout = {NOE,BINV,ncol};
    case 4,
        varargout = {NOE,BINV,ncol,NOBIN};
end
