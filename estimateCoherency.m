function [output,taper,concentration] = estimateCoherency(data1,data2,sampleFrequency,Ktapers,pad,freqsOfInterest,taper,concentration);
% function [output,taper,concentration] = estimateCoherency(data1,data2,sampleFrequency,Ktapers,pad,freqsOfInterest,taper,concentration);
%
% This function attempts to find the transfer function from data
[Nrepeats,Nsamples]=size(data1);
if nargin<3 sampleFrequency = []; end;
if nargin<4 Ktapers = []; end;
if nargin<5 pad = []; end;
if nargin<6 freqsOfInterest = []; end;
if nargin<7 taper = [];end;

if isempty(sampleFrequency) sampleFrequency = 1e3; end;
if isempty(Ktapers)  Ktapers = 2*ceil(Nsamples/sampleFrequency) - 1; end
if isempty(pad)  pad = max(2.^min(max((nextpow2(Nsamples)+4),12),16),Nsamples); end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFrequency/2]; 
elseif (length(freqsOfInterest)==1) freqsOfInterest = [0 freqsOfInterest]; end
bandwidth	= (Ktapers+1)*sampleFrequency/Nsamples; % The limit on frequency resolution is governed by the number of tapers
NW = (Ktapers+1)/2;
if isempty(taper); taper=dpss(Nsamples,NW); taper=taper(:,1:Ktapers); end;

freqsOfInterest = sort(freqsOfInterest);
pad_freqs=((1:pad)-1)*sampleFrequency/pad;
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));
freq_grid = pad_freqs(keep);
Nkeep = length(keep);

data1 = detrend(data1')';
data2 = detrend(data2')';

% Begin with eigenspectra
[xk1,sigPow1] = taperedSpectralEstimate(data1,taper,pad,1./sampleFrequency,keep,0); 
[xk2,sigPow2] = taperedSpectralEstimate(data2,taper,pad,1./sampleFrequency,keep,0);

% First look for lines in both processes
taperPower = sum(taper(:,1:2:Ktapers),1); % even tapers are nearly balanced, so only worry about
totalTaperPower = sum(taperPower.^2,2);   % odd ones for estimating power at each frequency

dims=size(xk1);
nd = length(dims);
% If we want to make maximal use of the independence of tapered estimates, we
% can take the last dimension, containing each of the tapered estimates, and
% bring it into the 1st dimension, the number of trials
xkt1 = reshape(permute(xk1,[nd 1:(nd-1)]),[dims(1)*dims(nd) dims(2:(nd-1))]);
for ch = 1:Nrepeats
    if (Ktapers>2)
        cs1(ch,:) = xkt1.'*taperPower'/totalTaperPower)'; % Least-squares estimate
        keyboard;
    else
        cs1(ch,:) = xk1(ch,:,1)*taperPower'/totalTaperPower;
    end
end
% Combine over trials and tapers to find line elements
totalPowerInBand = sum(squeeze(sum(abs(xk1.^2),3)),1);
explainedPower = sum(abs(cs1).^2,1);
residualPower = (totalPowerInBand/totalTaperPower-explainedPower);
fs1 = (Ktapers-1)*explainedPower./max(residualPower,eps);

totalPowerInBand = squeeze(sum(abs(xk2.^2),3));
for ch = 1:Nrepeats 
    if (Ktapers>2)
        cs1(ch,:) = (squeeze(xk2(ch,:,1:2:Ktapers))*taperPower'/totalTaperPower)'; % Least-squares estimate
    else
        cs1(ch,:) = xk2(ch,:,1)*taperPower'/totalTaperPower;
    end
end
explainedPower = abs(cs2).^2;
residualPower = (totalPowerInBand/totalTaperPower-explainedPower);
fs2 = (Ktapers-1)*explainedPower./max(residualPower,eps);

p1 = fcdf(fs1,2,2*(Ktapers-1));
p2 = fcdf(fs2,2,2*(Ktapers-1));

sig1 = p1>1-1/Bsamples;
sig2 = p2>1-1/Bsamples;
jointSig = find(sig1.*sig2);

