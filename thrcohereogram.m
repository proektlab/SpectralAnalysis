function [output,v,lambda] = thrcohereogram(TS1,TS2,sampleFreq,kTapers,freqsOfInterest,window,winstep,v,lambda)
% function [output,v,lambda] = thrcohereogram(TS1,sampleFreq,kTapers,freqsOfInterest,window,winstep,v,lambda)

if nargin<3 sampleFreq = []; end;
if nargin<4 kTapers = []; end;
if nargin<5 freqsOfInterest = []; end;
if nargin<6 window = []; end;
if nargin<7 winstep = []; end;
if nargin<8 v = []; end;
if nargin<9 lambda = [];end;

[nTrials,dataLength] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(kTapers)  kTapers = 5; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;
if isempty(window), window = nSamples; else nSamples = window; end;
if isempty(winstep), winstep = window/5; end;

sample_step	= 1/sampleFreq; 
fs			= sampleFreq/window; % Fundamental frequency in Hz
bandwidth	= (kTapers+1)*fs; % The limit on frequency resolution is governed by the number of tapers
timeResolution = window*sample_step / kTapers;
Tstep		= round(max(timeResolution / 2,0.02)*sampleFreq); %time window stepsize is 20 msec by default 
NW			=(window*sample_step)*(bandwidth/2); % Time-frequency product defines the limit on estimates
pad			=max(2^14,2^nextpow2(window)); % FFT runs fastest when a power of 2

if winstep < Tstep,
    Tstep = round(winstep);
end

midpoint = round(window/2);
ind = 0:Tstep:window-midpoint; % reconcile the window step and the time grid
winstep = ind(abs(ind - winstep)==min(abs(ind - winstep)));
ind=unique([midpoint:-Tstep:1 midpoint:Tstep:window]); % Create the time grid for the spectrograms
nTimePts=length(ind); % Number of time points in the high resolution grid

nSteps = round((dataLength-window)/winstep);
if (dataLength == window) 
    nSteps = 1; 
end;


[freqGrid,fkeep] = makeFrequencyGrid(sampleFreq,pad,freqsOfInterest,fs);

nFreqs = length(freqGrid);

if isempty(v)
   disp(' Calculating tapers');
   [v,lambda]=dpss(window,NW);			%generate the tapers
end
v=v(:,1:kTapers);
lambda=lambda(1:kTapers);


windowedData1 = zeros(nTrials*nSteps,window);
windowedData2 = windowedData1;
for j=1:nSteps
    inWindow = (1:window)+(j-1)*winstep;
    windowedData1((j-1)*nTrials+(1:nTrials),:) = TS1(:,inWindow);
    windowedData2((j-1)*nTrials+(1:nTrials),:) = TS2(:,inWindow);
end

dims = prod(size(windowedData1));
nRows = size(windowedData1,1);

brkpts = 1:ceil(1e7/dims*nRows):nRows;
if ~(brkpts(end) == size(windowedData1,1))
    brkpts = [brkpts size(windowedData1,1)];
elseif (length(brkpts)==1)
    brkpts(2) = 1;
end
for i = 1:length(brkpts)-1
    windowedData1(brkpts(i):brkpts(i+1),:) = detrend(windowedData1(brkpts(i):brkpts(i+1),:)')'; % Detrend after windowing the data
    windowedData2(brkpts(i):brkpts(i+1),:) = detrend(windowedData2(brkpts(i):brkpts(i+1),:)')'; % Detrend after windowing the data
end

w = waitbar(0,'Sliding window');
for i = 1:length(brkpts)-1
    waitbar(i/(length(brkpts)-1),w);
    nInBlock = 1+brkpts(i+1)-brkpts(i);
    [xk1,pow1] = taperedSpectralEstimate(windowedData1(brkpts(i):brkpts(i+1),:),v,pad,1/sampleFreq,fkeep,1);
    [xk2,pow2] = taperedSpectralEstimate(windowedData2(brkpts(i):brkpts(i+1),:),v,pad,1/sampleFreq,fkeep,1);

    if kTapers > 1
        d1 = computeAdaptiveWeights(xk1,pow1,pad,lambda);
        d1  = repmat(reshape(d1, [nInBlock,1,nFreqs,kTapers]),[1 nTimePts 1 1]);
        
        d2 = computeAdaptiveWeights(xk2,pow2,pad,lambda);
        d2  = repmat(reshape(d2, [nInBlock,1,nFreqs,kTapers]),[1 nTimePts 1 1]);
    else
        d1 = ones([nInBlock nTimePts nFreqs 1]);
        d2 = d1;
    end
    xk1 = repmat(reshape(xk1,[nInBlock,1,nFreqs,kTapers]),[1 nTimePts 1 1]);
    xk2 = repmat(reshape(xk2,[nInBlock,1,nFreqs,kTapers]),[1 nTimePts 1 1]);

    vGrid = repmat(reshape(v(ind,:),[1 nTimePts 1 kTapers]),[nInBlock 1 nFreqs 1]);

    % Back project on to time grid
    xk1 = xk1.*vGrid;
    xk2 = xk2.*vGrid;
    fisher1(brkpts(i):brkpts(i+1),:,:) = (mean(abs(xk1).^2,4)).^-2;
    fisher2(brkpts(i):brkpts(i+1),:,:) = (mean(abs(xk2).^2,4)).^-2;
    Sii(brkpts(i):brkpts(i+1),:,:) = computeSpectra(xk1,d1,xk1,d1,lambda);
    Sjj(brkpts(i):brkpts(i+1),:,:) = computeSpectra(xk2,d2,xk2,d2,lambda);
    Sij(brkpts(i):brkpts(i+1),:,:) = computeSpectra(xk1,d1,xk2,d2,lambda);
end
close(w);
indexTimes = ind - round(window/2);
baseTimes = round(window/2) + (0:nSteps-1)*winstep;
absoluteTimes = repmat(indexTimes,[nSteps 1])+repmat(baseTimes',[1 length(ind)]);
absoluteTimes = reshape(permute(absoluteTimes,[2 1]),[nSteps*nTimePts 1]);
timeGrid = unique(absoluteTimes);
K = reshape(permute(repmat(epanechnikov(nTimePts), [nRows 1]),[2 1]),[nSteps*nTimePts 1]);
avMat = zeros([length(timeGrid) nSteps*nTimePts]);
for i = 1:length(timeGrid), 
    avMat(i,:) = K.*(absoluteTimes == timeGrid(i)); 
end;
avMat = avMat./repmat(sum(avMat,2),[1 nSteps*nTimePts]);

tmp1 = reshape(permute(Sii,[2 1 3]),[nSteps*nTimePts nFreqs]);
spect1 = avMat*tmp1;

tmp2 = reshape(permute(Sjj,[2 1 3]),[nSteps*nTimePts nFreqs]);
spect2 = avMat*tmp2;
tmp12 = reshape(permute(Sij,[2 1 3]),[nSteps*nTimePts nFreqs]);
altcoh = tmp12./(sqrt(tmp1).*sqrt(tmp2));
spect12 = avMat*altcoh;

fisher1 = reshape(permute(fisher1,[2 1 3]),[nSteps*nTimePts nFreqs]);
fisher2 = reshape(permute(fisher2,[2 1 3]),[nSteps*nTimePts nFreqs]);

output.spect1 = spect1;
output.spect2 = spect2;
output.coh = spect12;
output.timeGrid = timeGrid/sampleFreq;
output.freqGrid = freqGrid;
output.bandwidth = bandwidth;
output.Ktapers = kTapers;
output.timeResolution = timeResolution;
output.winstep = winstep;
output.window = window;
output.fs = fs;
output.NW = NW;

end

function S = computeSpectra(xk1,d1,xk2,d2,lambda)
shape = size(xk1);
kTapers = length(lambda);
if shape(end) ~= kTapers
    shape = [shape 1];
end

A = sum(1./lambda);
lastDim = length(shape);

lambdaGrid = repmat(shiftdim(sqrt(lambda),1-lastDim),[shape(1:end-1) 1]);
S  = (A./kTapers)*sum(lambdaGrid.*d1.*xk1.*conj(lambdaGrid.*d2.*xk2),lastDim)./...
     (sqrt(sum(d1.^2,lastDim)).*sqrt(sum(d2.^2,lastDim)));
end

function [freqGrid,pad_fkeep2] = makeFrequencyGrid(sampleFreq,pad,freqsOfInterest,fs)
freqsOfInterest = sort(freqsOfInterest);
pad_freqs = sampleFreq*(1:pad)/pad;
pad_fkeep2 = find(pad_freqs>freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));

if fs<1 & freqsOfInterest(1)<.5 & freqsOfInterest(2)>80
    transitionpt = find(pad_freqs<10, 1, 'last' );
    %pad_fkeep2 = unique([round(logspace(log10(min(pad_fkeep2)),log10(transitionpt),150)) round(linspace(transitionpt,max(pad_fkeep2),150))]);
    pad_fkeep2 = unique(round(logspace(log10(min(pad_fkeep2)),log10(max(pad_fkeep2)),300)));
else
    frqStep = round(length(pad_fkeep2)/240);
    pad_fkeep2 = pad_fkeep2(1):frqStep:pad_fkeep2(end);
end;

freqGrid=pad_freqs(pad_fkeep2);
end


function [d,wk] = computeAdaptiveWeights(xk,power,pad,lambda)
[nTrials,nKeep,kTapers] = size(xk);
d = ones(nTrials,nKeep,kTapers);
wk = d;
for t = 1:nTrials
    Sk = squeeze(abs(xk(t,:,:)).^2);
    S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
    
    tol=.0005*power(t)/pad; % Set tolerance to match PMTM
    mser=power(t)*(1-lambda);
    
    temp=zeros(nKeep,1);
    S1=zeros(nKeep,1);
    while sum(abs(S0-S1)/nKeep)>tol
        d(t,:,:)=(S0*ones(1,kTapers))./(S0*lambda'+ones(nKeep,1)*mser'); % Percival & Walden p368.
        wk(t,:,:)=(squeeze(d(t,:,:)).^2).*(ones(nKeep,1)*lambda');
        S1=sum(squeeze(wk(t,:,:))'.*Sk')./ sum(squeeze(wk(t,:,:))'); % Percival & Walden p369.
        temp=S1';
        S1=S0;
        S0=temp;
    end
end
end

function K = epanechnikov(N)
K = .75*(1-((2*(1:N))/(N+1) - 1).^2);
K(K<0)=0;
end