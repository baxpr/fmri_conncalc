function B = fourier_filter_basis(N,tr,low,high)

%% Use dftmtx to generate a Fourier basis set for bandpass filtering
%
%    N      : number of points in the time series
%    low    : low freq cutoff in Hz
%    high   : high freq cutoff in Hz
%    tr     : sampling time (e.g. fMRI TR) in sec
%
%
% Even N:
%
%    real( D(:,1) )  is constant
%    imag( D(:,1) )  is constant
%    real( D(:,2) )  is lowest cos
%    imag( D(:,2) )  is lowest sin
%
%    real( D(:,floor(N/2)+1) )  is highest cos
%    imag( D(:,floor(N/2)+1) )  is constant
%
%    So we'd keep
%        cos: real(D(2:floor(N/2)+1))
%        sin: imag(D(2:floor(N/2)))
%
% Odd N:
%
%    real( D(:,1) )  is constant
%    imag( D(:,1) )  is constant
%    real( D(:,2) )  is lowest cos
%    imag( D(:,2) )  is lowest sin
%
%    real( D(:,floor(N/2)+1) )  is highest cos
%    imag( D(:,floor(N/2)+1) )  is highest sin
%
%    So we'd keep
%        cos: real(D(2:floor(N/2)+1))
%        sin: imag(D(2:floor(N/2)+1))
%
% The number of sin and cos terms in total should be N-1. The Nth term is
% the constant.

% Sampling rate in Hz
sr = 1 / tr;

% Frequency resolution
df = sr / N;

% Basis set
D = dftmtx(N);

% Select the non-redundant part
isodd = mod(N,2);
if isodd==0
	C = real(D(:,2:floor(N/2)+1));
	S = imag(D(:,2:floor(N/2)));
	fC = df * (1:floor(N/2));
	fS = df * (1:floor(N/2)-1);
else
	C = real(D(:,2:floor(N/2)+1));
	S = imag(D(:,2:floor(N/2)+1));
	fC = df * (1:floor(N/2));
	fS = df * (1:floor(N/2));
end

% Make bandpass matrix
bC = C(:,fC<low | fC>high);
bS = S(:,fS<low | fS>high);
B = [bC bS];

return


%% Diagnostics/plots

% Apply filter
X = [ones(size(B,1),1) B];
[~,~,fY] = regress(Y,X);

% Sample times
T = 0:tr:tr*(N-1);

% Plot
figure(1); clf; hold on
plot(T,Y,'-')
plot(T,fY,'-','LineWidth',2)

figure(2)
pwelch(fY,[],[],[],sr)
set(gca,'XTick',0:10:250)
