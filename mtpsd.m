function [output,taper,concentration]= mtpsd(TS1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration,doAdapt)
% function [output,taper,concentration]= mtpsd(TS1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration,doAdapt)
%
%
%
%
% 

if nargin<2 sampleFreq = []; end;
if nargin<3 Ktapers = []; end;
if nargin<4 freqsOfInterest = []; end;
if nargin<5 NW = []; end;
if nargin<6 pad = []; end;
if nargin<7 taper = []; end;
if nargin<8 concentration = []; end
if nargin<9 doAdapt = []; end

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(Ktapers)  Ktapers = 5; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;
if isempty(NW) NW = (Ktapers+1)/2; end;
if isempty(pad)  pad = 2.^max(nextpow2(nSamples)+4,12); end;
if isempty(taper); [taper,concentration]=dpss(nSamples,NW,'calc'); taper=taper(:,1:Ktapers); concentration=concentration(1:Ktapers); end;
if isempty(doAdapt) doAdapt = 1; end;



sample_step	= 1/sampleFreq; 
trials		= size(TS1,1); % Number of trials in the dataset
T				= size(TS1,2); % Length of a trial in sample points; determines analysis window
bandwidth	= (Ktapers+1)/(T*sample_step); % The limit on frequency resolution is governed by the number of tapers
fs				=1/(T*sample_step); % Fundamental frequency in Hz
total_harms	=(1/(2*sample_step))/fs; % Nyquist divided by frequency of 1st harmonic (fundamental frequency)
all_freqs	=fs*(1:total_harms); % Define frequency grid
hgrid			=1:total_harms; % Points in frequency grid
freqs_on_hgrid	=fs:bandwidth:(total_harms*fs); % Sampling determined by frequency grid
NW				=(T*sample_step)*(bandwidth/2); % Time-frequency product defines the limit on estimates

freqsOfInterest = sort(freqsOfInterest);
pad_freqs = sampleFreq*(0:1/pad:1-1/pad);
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));
 if pad ~= T
     frqStep = ceil(length(keep)/500);
     keep = keep(1):frqStep:keep(end);
 end
output.freq_grid = pad_freqs(keep);

TS1=detrend(TS1')';

rand('seed',sum(clock));

disp('Begining Power Calculation');
[dft,power] = taperedSpectralEstimate(TS1,taper,pad,1/sampleFreq,keep,0);
dft = suppressEigencoefficientSidelobes(dft,power,taper,concentration,pad);
psd = mean(2.*dft.*conj(dft)/sampleFreq,3);
%%% Compute confidence intervals
alpha = 0.05;
dof = 2*Ktapers*trials; % Thomson and Chave [in Advances in Spectrum Analysis and Array Processing, vol 1], p76;
q = chi2inv([alpha/2 1-alpha/2],dof);

if trials==1 &~doAdapt
    [CI,SE] = jackSpectralCI(permute(squeeze(psd),[2 1]),alpha);
elseif trials==1 & doAdapt
    [CI,SE] = jackSpectralCI(permute(squeeze(2*dft.*conj(dft)/sampleFreq),[2 1]),alpha);
elseif size(psd,3)>1
    [CI,SE] = jackSpectralCI(mean(psd,3),alpha);
else
    [CI,SE] = jackSpectralCI(psd,alpha);
end
output.xk=dft;
output.psd=psd;
output.Ktapers = Ktapers;
output.bandwidth = bandwidth;
output.jackknifeCI = CI;
output.jackStdErrdB = SE;
output.chi2CI = [dof*mean(psd,1)/q(2); dof*mean(psd,1)/q(1)];
return;

%%%%%%%%%%%
% Helper function for computing power spectral density confidence limits;
%
% this is based on using the log transform of the jackknifed spectral 
% estimates to develop nearly normal random variables, and then defining
% the confidence limits on the jackknife estimate of the variance of this
% quantity.  See Thomson and Chave, 1991
%%%%%%%%%%%
function [CI,jkStdErr] = jackSpectralCI(data,alpha);
if (nargin<2)
    alpha = 0.05;
end
N = size(data,1);
CInum=tinv(1-alpha/2,N-1);
jkLogSpecHat = log((1-eye(N))*data/(N-1));
jkMean = mean(jkLogSpecHat,1);
jkStdErr = sqrt(((N-1)/N)*sum((jkLogSpecHat-repmat(jkMean,[N 1])).^2,1));
CI = [mean(data,1).*exp(jkStdErr.*CInum); mean(data,1).*exp(-jkStdErr.*CInum)];
jkStdErr = 10*log10(exp(1))*jkStdErr; % rescale to dBs
return;
