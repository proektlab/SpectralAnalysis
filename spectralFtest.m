function [fs,cs,freq_grid,taper] = spectralFtest(data,sample_step,Ktapers,pad,freqsOfInterest,taper)
% function [fs,cs,freq_grid] = spectralFtest(data,sample_step,Ktapers,pad,freqsOfInterest,taper)
%       fspectrum, the fstatistics of the spectrum.  Asymptotically
%           distributed as F_2_,_2(K-1).
%       cs, the complex spectrum; can be used to reconstruct the line
%           elements
%       freq_grid, a vector of the frequencies contained in the spectral
%           quantities
%
% Written by Andrew Hudson, 7/16/2006. Version 1.3, last edited 6/23/08.

[Nchannels,Nsamples]=size(data);
if nargin<2, sample_step = []; end;
if nargin<3, Ktapers = []; end;
if nargin<4, pad = []; end;
if nargin<5, freqsOfInterest = []; end;
if nargin<6, taper = [];end;

if isempty(sample_step), sample_step = 0.001; end; sampleFreq = 1./sample_step;
if isempty(Ktapers),  Ktapers = 2*ceil(Nsamples/sampleFreq) - 1; end
if isempty(pad),  pad = max(2.^min(max((nextpow2(Nsamples)+4),12),16),Nsamples); end
if isempty(freqsOfInterest), freqsOfInterest = [0 1/(2*sample_step)]; 
elseif (length(freqsOfInterest)==1), freqsOfInterest = [0 freqsOfInterest]; end
NW = Ktapers+1;
if isempty(taper), [taper]=dpss(Nsamples,NW); taper=taper(:,1:Ktapers); end;

freqsOfInterest = sort(freqsOfInterest);
pad_freqs=(1:pad)/(pad*sample_step);
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));
freq_grid = pad_freqs(keep);
data = detrend(data')';
xk = taperedSpectralEstimate(data,taper,pad,sample_step,keep,0);

taperPower = sum(taper(:,1:2:Ktapers),1); % even tapers are nearly balanced, so we only worry about odd ones
totalTaperPower = sum(taperPower.^2,2);
bandPower = squeeze(sum(abs(xk.^2),3));
cs = zeros(Nchannels,size(xk,2));
for ch = 1:Nchannels % this could be vectorized, but for memory's sake we'll leave it be
    if (Ktapers>2)
        temp = squeeze(xk(ch,:,1:2:Ktapers));
        temp = temp*taperPower';
        cs(ch,:) = temp/totalTaperPower';
    else
        temp = xk(ch,:,1)*taperPower';
        cs(ch,:) = temp/totalTaperPower;
    end
end
linePower = abs(cs).^2;
fs = (Ktapers-1)*linePower./max((bandPower/totalTaperPower-linePower),eps);