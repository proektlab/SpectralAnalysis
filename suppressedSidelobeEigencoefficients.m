function [output,taper,concentration]= suppressedSidelobeEigencoefficients(TS1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration)
% function [xk,taper,concentration]= suppressedSidelobeEigencoefficients(TS1,sampleFreq,Ktapers,freqsOfInterest,NW,pad,taper,concentration)
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
output.freq_grid = pad_freqs(keep);

TS1=detrend(TS1')';

rand('seed',sum(clock));

disp('Begining Power Calculation');
[dft,power] = taperedSpectralEstimate(TS1,taper,pad,1/sampleFreq,keep,0);

if (Ktapers>1) & ~isempty(concentration)
    w2 = waitbar(0,'Adaptive weighting of estimates');
    nKeep = length(keep);
    for t = 1:size(dft,1)
        waitbar(t/(size(dft,1)),w2);
        Sk = squeeze(abs(dft(t,:,:)).^2);
        S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
        
        tol=.0005*power(t); % Set tolerance to match PMTM
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
        output.xk(t,:,:) = (sqrt(ones(nKeep,1)*concentration').*repmat(S0,[1 Ktapers]))...
            ./(ones(nKeep,1)*concentration'.*repmat(S0,[1 Ktapers]) + (ones(nKeep,1)*mser')).*squeeze(dft(t,:,:));
        output.psd(t,:) = 2*S0'/sampleFreq;
    end
    close(w2);
else
    disp('Too few tapers to do adaptive fitting; returning raw estimates.');
    output.xk = dft;
    output.psd = 2*dft.*conj(dft)/sampleFreq;
end
return;
