function [output,v,lambda]= thrautocoherence(TS1,sampleFreq,kTapers,freqsOfInterest,v,lambda)

if nargin<2 sampleFreq = []; end;
if nargin<3 kTapers = []; end;
if nargin<4 freqsOfInterest = []; end;
if nargin<5 v = []; end;
if nargin<6 lambda = [];end;

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(kTapers)  kTapers = 5; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;

sample_step	= 1/sampleFreq; 
T			= size(TS1,2); % Length of a trial in sample points; determines analysis window
fs			= sampleFreq/T; % Fundamental frequency in Hz
bandwidth	= (kTapers+1)*fs; % The limit on frequency resolution is governed by the number of tapers
NW			=(T*sample_step)*(bandwidth/2); % Time-frequency product defines the limit on estimates
pad			=max(2^14,2^nextpow2(T)); % FFT runs fastest when a power of 2


[use_freqs,fkeep] = makeFrequencyGrid(sampleFreq,pad,freqsOfInterest,fs);

nFreqs = length(use_freqs);

if isempty(v)
   disp(' Calculating tapers');
   [v,lambda]=dpss(T,NW);			%generate the tapers
end
v=v(:,1:kTapers);
lambda=lambda(1:kTapers);

TS1 = detrend(TS1')';

[xk1,pow1] = taperedSpectralEstimate(TS1,v,pad,sample_step,fkeep,1);
d1 = computeAdaptiveWeights(xk1,pow1,pad,lambda);

xk1 = repmat(reshape(xk1,[nTrials,1,nFreqs,kTapers]),[1 nFreqs 1 1]);
d1  = repmat(reshape(d1, [nTrials,1,nFreqs,kTapers]),[1 nFreqs 1 1]);

for f1 = 1:nFreqs
    Sij(:,:,f1) = computeSpectra(reshape(xk1(:,f1,:,:),[nTrials nFreqs kTapers]),...
        reshape(d1(:,f1,:,:),[nTrials nFreqs kTapers]),...
        reshape(xk1(:,:,f1,:),[nTrials nFreqs kTapers]),...
        reshape(d1(:,:,f1,:),[nTrials nFreqs kTapers]),lambda);
end
Sii = zeros(nTrials,nFreqs);
SiiSjj = zeros(nTrials,nFreqs,nFreqs);
for t=1:nTrials,
    Sii(t,:) = diag(squeeze(Sij(t,:,:)));
    SiiSjj(t,:,:) = sqrt(Sii(t,:)')*sqrt(Sii(t,:));
    SiiSjj(t,:,:) = (1-eye(nFreqs)).*squeeze(SiiSjj(t,:,:)) + diag(Sii(t,:));
end

output.autocoh = Sij./SiiSjj;
output.spect = Sii;
output.freqGrid = use_freqs;
output.bandwidth = bandwidth;

end

function [S,dof] = computeSpectra(xk1,d1,xk2,d2,lambda);
shape = size(xk1);
lastDim = length(shape);
kTapers = length(lambda);
if shape(lastDim) ~= kTapers
    shape = [shape 1];
    lastDim = length(shape);
end

lambdaGrid = repmat(shiftdim(sqrt(lambda),1-lastDim),[shape(1:end-1) 1]);

A = sum(1./lambda);
S = (A./kTapers)*sum(lambdaGrid.*d1.*xk1.*conj(lambdaGrid.*d2.*xk2),lastDim)./...
    (sqrt(sum(d1.^2,lastDim)).*sqrt(sum(d2.^2,lastDim)));

end

function [freqGrid,pad_fkeep2] = makeFrequencyGrid(sampleFreq,pad,freqsOfInterest,fs)
freqsOfInterest = sort(freqsOfInterest);
pad_freqs = sampleFreq*(1:pad)/pad;
pad_fkeep2 = find(pad_freqs>freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));

%if fs<1 & freqsOfInterest(1)<.5 & freqsOfInterest(2)>80
if fs<1 & freqsOfInterest(1)<.5
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