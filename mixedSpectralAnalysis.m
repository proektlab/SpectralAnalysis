function [output,taper,concentration]= mixedSpectrumAnalysis(data,F1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration)
% function [xk,taper,concentration]= suppressedSidelobeEigencoefficients(TS1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration)
%
%
%
%
% 

% Preliminaries
if nargin<2 sampleFreq = []; end;
if nargin<3 Ktapers = []; end;
if nargin<4 freqsOfInterest = []; end;
if nargin<5 NW = []; end;
if nargin<6 pad = []; end;
if nargin<7 taper = []; end;
if nargin<8 concentration = []; end

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) 
    sampleFreq = 1e3; 
end;
if isempty(Ktapers)  
    Ktapers = 5; 
end
if isempty(freqsOfInterest) 
    freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 
    freqsOfInterest = [0 freqsOfInterest]; 
end;
if isempty(NW) 
    NW = (Ktapers+1)/2; 
elseif Ktapers > 2*NW -1
    warning('You have requested more tapers than is possible at this time-frequency resolution.');
end;
if isempty(pad)  
    pad = 2.^max(nextpow2(nSamples)+6,17);
end;
if isempty(taper); 
    [taper,concentration]=dpss(nSamples,NW,'calc'); 
    taper=taper(:,1:Ktapers); 
    concentration=concentration(1:Ktapers); 
end;
% done with preliminaries

sampleStep	= 1/sampleFreq; 
[Ntrials,Nsamples] = size(TS1);
frequencyBandwidth = (Ktapers+1)*sampleFreq/Nsamples; % this is full-width, 2W
rayleighFrequency = sampleFreq/Nsamples;
if F1<frequencyBandwidth,
    warning('The fundamental is not well resolved from DC at this frequency bandwidth.  The regression for line elements may become ill-conditioned.');
end;
if F1<rayleighFrequency
    warning('The fundamental is less than the Rayleigh frequency. Are you sure you want to do this?');
end


freqsOfInterest = sort(freqsOfInterest);
meshGrid = sampleFreq*(0:1/pad:1-1/pad);
keep = find(meshGrid>=freqsOfInterest(1)&meshGrid<=freqsOfInterest(2));
output.freq_grid = meshGrid(keep);


DC = mean(swatch,2);
% Generate Eigenspectra
[xjk,sigmaj]=taperedSpectralEstimate(data,taper,pad,1/sampleFreq,keep,1);

