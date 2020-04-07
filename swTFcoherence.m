function [out,taper] = swTFcoherence(TS1,TS2,sampleFreq,Ktapers,freqsOfInterest,window,winstep,NW,pad,taper);
% function [out,taper] = swTFcoherence(TS1,TS2,sampleFreq,Ktapers,freqsOfInterest,window,winstep,NW,pad,taper);
%
% INPUT:
%       TS1, the [nTrials, nSamples] matrix of data from time series 1
%       TS2, the [nTrials, nSamples] matrix of data from time series 2
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
% OUTPUT:
%       output, a structure with the following fields:
%           tfse1, the time-frequency spectrogram for the first timeseries
%           tfse2, the time-frequency spectrogram for the second timeseries
%           tfcoh, the time-frequency cohereogram for the two timeseries
%           freq_grid, the frequencies in the spectrogram
%           time_grid, the timepoints in the spectrogram
%           bandwidth, the bandwith of the spectrogram
%
if nargin<3 sampleFreq = []; end;
if nargin<4 Ktapers = []; end;
if nargin<5 freqsOfInterest = []; end;
if nargin<6 window = []; end;
if nargin<7 winstep = [];end;
if nargin<8 NW = []; end;
if nargin<9 pad = []; end;
if nargin<10 taper = [];end;

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(Ktapers)  Ktapers = 2*ceil(nSamples/sampleFreq) - 1; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;
if isempty(window) window = min(250,nSamples); end;
if isempty(winstep) winstep = window/5; end;
if isempty(NW) NW = Ktapers+1; end;
if isempty(pad)  pad = 2.^max(nextpow2(nSamples+4),12); end;
if isempty(taper); taper=dpss(window,NW,'calc'); taper=taper(:,1:Ktapers); end;

baseFreq = 1./(window/sampleFreq);

freqsOfInterest = sort(freqsOfInterest);
pad_freqs = sampleFreq*(1:pad)/pad;
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));

if baseFreq<.5 && freqsOfInterest(1)<.5 && freqsOfInterest(2)>80
    transitionpt = find(pad_freqs<10, 1, 'last' );
    keep = unique([round(logspace(log10(min(keep)),log10(transitionpt),150)) round(linspace(transitionpt,max(keep),150))]);
else
    frqStep = round(length(keep)/240);
    keep = keep(1):frqStep:keep(end);
end;

nSteps = round((nSamples-window)/winstep);
taper = taper(:,1:Ktapers);

out.time_grid = linspace(0,nSamples/sampleFreq,nSteps);
out.freq_grid = sampleFreq*(1:pad/2)./pad;
out.freq_grid = out.freq_grid(keep);


windowedTS1 = zeros(nTrials*nSteps,window);
windowedTS2 = zeros(nTrials*nSteps,window);
for j=1:nSteps
    inWindow = (1:window)+(j-1)*winstep;
    windowedTS1((j-1)*nTrials+(1:nTrials),:) = TS1(:,inWindow);
    windowedTS2((j-1)*nTrials+(1:nTrials),:) = TS2(:,inWindow);
end
windowedTS1 = detrend(windowedTS1')';
windowedTS2 = detrend(windowedTS2')';

brkpts = [1:10000:size(windowedTS1,1)];
if ~(brkpts(end) == size(windowedTS1,1))
     brkpts = [brkpts size(windowedTS1,1)];
 end

tfse1 = zeros(brkpts(end),length(keep));
tfse2 = tfse1;
tfcoh = tfse1;

w = waitbar(0,'Sliding window');
for i = 1:length(brkpts)-1
    waitbar(i/(length(brkpts)-1),w);
    s1 = taperedSpectralEstimate(windowedTS1(brkpts(i):brkpts(i+1),:),taper,pad,1/sampleFreq,keep,1);
    s2 = taperedSpectralEstimate(windowedTS2(brkpts(i):brkpts(i+1),:),taper,pad,1/sampleFreq,keep,1);
    temp1 = sum(s1.*conj(s1),3);
    tfse1(brkpts(i):brkpts(i+1),:) = temp1/Ktapers;
    temp2 = sum(s2.*conj(s2),3);
    tfse2(brkpts(i):brkpts(i+1),:) = temp2/Ktapers;
    temp3 = sum(s1.*conj(s2),3);
    tfcoh(brkpts(i):brkpts(i+1),:) = temp3./sqrt(abs(temp1).*abs(temp2));;
end
close(w);

for j=1:nSteps
    out.tfse1(:,j,:)  = tfse1((j-1)*nTrials+(1:nTrials),:);
    out.tfse2(:,j,:)  = tfse2((j-1)*nTrials+(1:nTrials),:);
    out.tfcoh(:,j,:)  = tfcoh((j-1)*nTrials+(1:nTrials),:);
end
clear tfse;
out.bandwidth = (Ktapers+1)*sampleFreq/window;

[whereSigUP,sigCohUP,uniformPhaseP]=uniformPhaseTest(out.tfcoh,0.01);
out.uniformPhaseP = uniformPhaseP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [magCI,phaseCI,magP,phaseP] = jackCoherenceCI(s11,s22,s12,criterion);
% This function uses the variance stabilizing transform of Thomson and Chave, 1991.
if nargin<4
   criterion=0.05
end
nd = ndims(s11);
dims = size(s11);
if nd>2
   s11 = reshape(permute(s11,[nd 1:(nd-1)]),[dims(1)*dims(nd) dims(2:(nd-1))]);
   s22 = reshape(permute(s22,[nd 1:(nd-1)]),[dims(1)*dims(nd) dims(2:(nd-1))]);
   s12 = reshape(permute(s12,[nd 1:(nd-1)]),[dims(1)*dims(nd) dims(2:(nd-1))]);
   dims=size(s11);
end
N = dims(1);
CInum = tinv(1-criterion/2,N-1);
jkMat = 1-eye(N);
jkC = jkMat*s12./sqrt((jkMat*s11).*(jkMat*s22));
jkQ = sqrt(2*N-2)*atanh(abs(jkC));
jkMeanQ = mean(jkQ,1);
jkStdErrQ = sqrt(((N-1)/N)*sum((jkQ-repmat(jkMeanQ,[N 1])).^2,1));
magCI = [jkMeanQ+CInum.*jkStdErrQ; jkMeanQ-CInum.*jkStdErrQ];
magCI = tanh(magCI./sqrt(2*N-2));
magP = normcdf(0,jkMeanQ,jkStdErrQ);

jkE = jkC./abs(jkC);   % phase factor
jkMeanE = mean(jkE,1); 
jkStdErrPhi = sqrt(2*(N-1)*(1-abs(jkMeanE))); % Thomson and Chave (1991), p91, eq 2.63
%meanPhi = angle(mean(s11,1)./sqrt(mean(s11,1).*mean(s22,1)));
%phaseCI = [meanPhi + CInum.*jkStdErrPhi; meanPhi - CInum.*jkStdErrPhi];
phaseCI = [CInum.*jkStdErrPhi; - CInum.*jkStdErrPhi];
phaseP = normcdf(0,jkMeanE,jkStdErrPhi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [muHat,sigmaHat,jk] = jackknife(data,jackOverTapers)
% This function computes the mean and standard deviation of coherency in 
% a dataset by using a jackknife.
% Each row of the data is assumed to be an independent trial.  If the 
% logical flag, jackOverTapers, is set, the jackknife will act along 
% the first and last dimensions of the data.

nd = ndims(data);

if nargin<2 | nd<3
   jackOverTapers=0;
end

dims=size(data);

% If we want to make maximal use of the independence of tapered estimates, we
% can take the last dimension, containing each of the tapered estimates, and
% bring it into the 1st dimension, the number of trials
if jackOverTapers
   data = reshape(permute(data,[nd 1:(nd-1)]),[dims(1)*dims(nd) dims(2:(nd-1))]);
   dims = size(data); % We've changed the size of the data
end
N = dims(1);
% Create jackknife population
jk = zeros(N,prod(dims(2:end)));
if length(dims)>2
   jk = reshape((1-eye(N))*reshape(data, [N prod(dims(2:end))])/(N-1),dims);
else
   jk = (1-eye(N))*data./(N-1);
end

% Make relevant calculations
muHat=mean(data,1);
sigmaHat=sqrt(((N-1)/N).*sum((jk - repmat(muHat,[N 1 1])).^2,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [muHat,sigmaHat] = jkVar(data)
% This function computes the mean and standard deviation of coherency in 
% a dataset by using a jackknife.
% Each row of the data is assumed to be an independent trial.  If the 
% logical flag, jackOverTapers, is set, the jackknife will act along 
% the first and last dimensions of the data.
N = size(data,1);
muHat = atanh(abs(mean(data,1)));
sigmaHat=sqrt(((N-1)/N).*sum((atanh(abs(data)) - repmat(muHat,[N 1 1])).^2,1));
