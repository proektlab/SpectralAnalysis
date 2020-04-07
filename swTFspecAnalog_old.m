function [out,taper,concentration] = swTFspecAnalog(data,sampleFreq,Ktapers,freqsOfInterest,window,winstep,NW,pad,taper,concentration,doAdapt,includePhase)
%function [out,taper,concentration] = swTFspecAnalog(data,sampleFreq,Ktapers,freqsOfInterest,window,winstep,NW,pad,taper,concentration,doAdapt)
%
% INPUT:
%       data, the [nTrials, nSamples] matrix of data
%       sampleFreq, the sample frequency in Hz
%       Ktapers, the number of tapers to use in the calculation 
%           (default: 2*ceil(Nsamples/sampleFreq)-1 )
%       freqsOfInterest, a two element vector of frequencies to return,
%           from min(freqsOfInterest) to max(freqsOfInterest)
%       window, the sliding window length in sample points (default: 250)
%       winstep, the number of data points to step between windows
%           (default: window/5)
%       NW, the time-frequency product (default: Ktapers + 1)
%       pad, the pad length for the fft (default: next power of 2>16*window)
%       taper, the taper to use (default: calculates the appropriate
%           Slepian tapers for the requested parameters)
%       concentration, the eigenvalues corresponding to the Slepian tapers;
%           this is required for adaptive side-lobe supresssion and not
%           used otherwise; if no tapers are provided, this will be
%           automatically calculated during the generation of the tapers
%       doAdapt, a logical flag to control whether to engage in adaptive
%           sidelobe supresssion to minimize estimation bias when using
%           Slepian tapers.  Defaults to true.
%
% OUTPUT:
%       output, a structure with the following fields:
%       tfse, the time-frequency spectrogram
%       freq_grid, the frequencies in the spectrogram
%       time_grid, the timepoints in the spectrogram
%       bandwidth, the bandwith of the spectrogram
%
% Written by Andrew Hudson, 10/17/2004.  Version 2.2 last updated 3/15/12.


if nargin<2, sampleFreq = []; end;
if nargin<3, Ktapers = []; end;
if nargin<4, freqsOfInterest = []; end;
if nargin<5, window = []; end;
if nargin<6, winstep = [];end;
if nargin<7, NW = []; end;
if nargin<8, pad = []; end;
if nargin<9, taper = [];end;
if nargin<10, concentration = []; end;
if nargin<11, doAdapt = []; end;
if nargin<12, includePhase = []; end;

[nTrials,nSamples] = size(data);
if nSamples == 1 && nTrials>1
    data = data';
    [nTrials,nSamples] = size(data);
end;    

if isempty(sampleFreq), sampleFreq = 1000; end;
if isempty(Ktapers),  Ktapers = 2*ceil(nSamples/sampleFreq) - 1; end  % try for 2 Hz bandwidth
if isempty(freqsOfInterest), freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1, freqsOfInterest = [0 freqsOfInterest]; end;
if isempty(window), window = min(250,nSamples); end;
if isempty(winstep), winstep = window/5; end;
if isempty(NW), NW = (Ktapers+1)/2; end;
if isempty(pad),  pad = 2.^max(nextpow2(window)+4,12); end; % 
if Ktapers>2*NW, warning(1,['Using the maximum of ' num2str(2*NW) ' tapers for the availible data.']); Ktapers = 2*NW; end;
if isempty(taper), [taper,concentration]=dpss(window,NW,'calc'); taper=taper(:,1:Ktapers); concentration=concentration(1:Ktapers); end;
if isempty(doAdapt), doAdapt = 1; end;
if isempty(includePhase), includePhase = 1; end;

freqsOfInterest = sort(freqsOfInterest);
baseFreq = 1./(window/sampleFreq);

pad_freqs = sampleFreq*(0:1/pad:1-1/pad);
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));

if baseFreq<.5 && freqsOfInterest(1)<.5 && freqsOfInterest(2)>80
    transitionpt = find(pad_freqs<10, 1, 'last' );
    keep = unique([round(logspace(log10(min(keep)),log10(transitionpt),150)) round(linspace(transitionpt,max(keep),150))]);
else
    frqStep = round(length(keep)/240);
    keep = keep(1):frqStep:keep(end);
end;

nSteps = round((nSamples-window)/winstep);
if (nSamples == window) 
    nSteps = 1; 
end;
taper = taper(:,1:Ktapers);

%out.time_grid = linspace(0,nSamples/sampleFreq,nSteps);
%out.time_grid = (window/2 + (0:winstep:nSamples-window-1))/sampleFreq;
out.time_grid = (window/2 + ((1:nSteps)-1)*winstep)/sampleFreq;
out.freq_grid = pad_freqs;
out.freq_grid = out.freq_grid(keep);


windowedData = zeros(nTrials*nSteps,window);
for j=1:nSteps
    inWindow = (1:window)+(j-1)*winstep;
    windowedData((j-1)*nTrials+(1:nTrials),:) = data(:,inWindow);
end
windowedData = detrend(windowedData')'; % Detrend after windowing the data

brkpts = 1:10000:size(windowedData,1);
if ~(brkpts(end) == size(windowedData,1))
    brkpts = [brkpts size(windowedData,1)];
elseif (length(brkpts)==1)
    brkpts(2) = 1;
end

tfse = zeros(brkpts(end),length(keep));
%dof = tfse;
if includePhase, phase = tfse; end;

w = waitbar(0,'Sliding window');
for i = 1:length(brkpts)-1
    waitbar(i/(length(brkpts)-1),w);
    [dft,power] = taperedSpectralEstimate(windowedData(brkpts(i):brkpts(i+1),:),taper,pad,1/sampleFreq,keep,0);
    if ((doAdapt && Ktapers>1) && ~isempty(concentration))
        w2 = waitbar(0,'Adaptive weighting of estimates');
        nKeep = length(keep);
        for t = 1:size(dft,1)
            if mod(t,10) == 0, waitbar(t/(size(dft,1)),w2); end;
            Sk = squeeze(abs(dft(t,:,:)).^2);
            S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
            
            tol=.0005*power(t)/pad; % Set tolerance to match PMTM
            mser=power(t)*(1-concentration);

            temp=zeros(nKeep,1);
            S1=zeros(nKeep,1);
            while sum(abs(S0-S1)/nKeep)>tol
                b=(S0*ones(1,Ktapers))./(S0*concentration'+ones(nKeep,1)*mser'); % Percival & Walden p368.
                wk=(b.^2).*(ones(nKeep,1)*concentration'); 
                S1=sum(wk'.*Sk')./ sum(wk'); % Percival & Walden p369.
                temp=S1'; 
                S1=S0; 
                S0=temp; 
            end
            tfse(brkpts(i)+t-1,:) = 2*S0'/sampleFreq;
            %dof(brkpts(i)+t-1,:) = 2*sum(wk.^2,2)./sum((b.^4).*(ones(nKeep,1)*concentration'),2);

        end
        close(w2);
    else
        temp = sum(dft.*conj(dft),3);
        tfse(brkpts(i):brkpts(i+1),:) = 2*temp/(Ktapers*sampleFreq);
    end
    if includePhase, phase(brkpts(i):brkpts(i+1),:) = angle(mean(dft,3)); end;
end
close(w);
clear temp dft power;
for j=1:nSteps
    out.tfse(:,j,:) = tfse((j-1)*nTrials+(1:nTrials),:);
end
if includePhase, 
    for j=1:nSteps
        out.phase(:,j,:) = phase((j-1)*nTrials+(1:nTrials),:); 
    end;
end

clear tfse phase;
out.bandwidth = (Ktapers+1)*sampleFreq/window/2;
out.Ktapers = Ktapers;
