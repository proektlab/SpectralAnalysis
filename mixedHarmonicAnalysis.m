function [out,taper,concentration] = mixedHarmonicAnalysis(data,sample_step,f1,Ktapers,pad,freqsOfInterest,taper,concentration)
% function [out,taper,concentration] = mixedHarmonicAnalysis(data,sample_step,f1,Ktapers,pad,freqsOfInterest,taper,concentration)

[Nchannels,Nsamples]=size(data);

if nargin<2 sample_step = []; end;
if nargin<3 f1 = []; end;
if nargin<4 Ktapers = []; end;
if nargin<5 pad = []; end;
if nargin<6 freqsOfInterest = []; end;
if nargin<7 taper = []; end;
if nargin<8 concentration = []; end;

if isempty(sample_step) sample_step = 0.001; end;
sampleFreq = 1/sample_step;
if isempty(Ktapers)  Ktapers = 2*ceil(Nsamples*sample_step) - 1; end
if isempty(pad)  pad = max(2.^min(max((nextpow2(Nsamples)+4),12),16),Nsamples); end
if isempty(freqsOfInterest) freqsOfInterest = [0 1/(2*sample_step)]; 
elseif (length(freqsOfInterest)==1) freqsOfInterest = [0 freqsOfInterest]; end
alpha = 0.05;
NW = Ktapers+1;

bandwidth	= NW/(sample_step*Nsamples); % The limit on frequency 
% resolution is governed by the number of tapers; this is
% Percival and Walden's 2W parameter.


if isempty(taper)|isempty(concentration),
    [taper,concentration]=dpss(Nsamples,NW,'calc'); 
    taper=taper(:,1:Ktapers); 
    concentration=concentration(1:Ktapers); 
end;

freqsOfInterest = sort(freqsOfInterest);
pad_freqs=(0:pad-1)/(pad*sample_step);
keep = find(pad_freqs>freqsOfInterest(1)&pad_freqs<freqsOfInterest(2));
freq_grid = pad_freqs(keep);

% Clear low frequencies
a = [(1:Nsamples)/Nsamples; ones(1,Nsamples)]';
linfit = (a*(a\data'))'; % the least-squares linear fit
data = data - linfit; 
DCoffset = mean(linfit,2);  % save the DC and linear components
linslope = (linfit(:,Nsamples) - linfit(:,1))./(Nsamples*sample_step);


[jk,sigPow] = taperedSpectralEstimate(data,taper,pad,sample_step,keep,0); % tapered spectral estimates arranged as [channel,frequency,taper]
H = fft(taper,pad);
taperPower = H(1,1:2:Ktapers); % even tapers are nearly balanced, so we only worry about odd ones
totalTaperPower = sum(taperPower.^2,2);
sp = squeeze(sum(abs(jk.^2),3));
for ch = 1:Nchannels % this could be vectorized, but for memory's sake we'll leave it be
    if (Ktapers>2)
        cs(ch,:) = (squeeze(jk(ch,:,1:2:Ktapers))*taperPower'/totalTaperPower)';
    else
        cs(ch,:) = jk(ch,:,1)*taperPower'/totalTaperPower;
    end
end
num = abs(cs).^2;
fs = (Ktapers-1)*num./max((sp/totalTaperPower-num),eps);

% Translate into psd units
jk = jk.*sqrt(sample_step);

% Search for significant components at harmonics of f1
harmonics = f1*(1:floor(freqsOfInterest(2)/f1));
harmonics = harmonics(find(harmonics>freqsOfInterest(1)));
Nharmonics = length(harmonics);
criterion = max(1-1/Nharmonics,0.95);

for ch = 1:Nchannels,
    Fstat(ch,:) = interp1(freq_grid,fs(ch,:),harmonics);
    Cstat(ch,:) = interp1(freq_grid,cs(ch,:),harmonics);
    Pval(ch,:) = fcdf(Fstat(ch,:),2,2*(Ktapers-1));
end
sig = Pval>criterion;

% Construct the mixed spectrum estimate
nKeep = length(keep);
if Ktapers>1
    w2 = waitbar(0,'Adaptive weighting of estimates for mixed spectrum');
    for t = 1:size(jk,1)
        waitbar(t/(size(jk,1)),w2);
        Sk = squeeze(abs(jk(t,:,:)).^2);
        S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
        
        tol=.0005*sigPow(t)/pad; % Set tolerance to match PMTM
        mser=sigPow(t)*(1-concentration);
        
        temp=zeros(nKeep,1);
        S1=zeros(nKeep,1);
        b=(S0*ones(1,Ktapers))./(S0*concentration'+ones(nKeep,1)*mser'); % Percival & Walden p368.
        wk=(b.^2).*(ones(nKeep,1)*concentration'); 
        while sum(abs(S0-S1)/nKeep)>tol
            b=(S0*ones(1,Ktapers))./(S0*concentration'+ones(nKeep,1)*mser'); % Percival & Walden p368.
            wk=(b.^2).*(ones(nKeep,1)*concentration'); 
            S1=sum(wk'.*Sk')./ sum(wk'); % Percival & Walden p369.
            temp=S1'; 
            S1=S0; 
            S0=temp; 
        end
        dof1(t,:) = (2*(sum(wk,2)).^2./sum(wk.^2,2))';
        mixedpsd(t,:) = 2*S0';
    end
    waitbar(1,w2,'Removing line elements');
else
    mixedpsd = abs(jk).^2;
    dof1 = 2;
end
q1 = [shiftdim(chi2inv([alpha/2],dof1),-1);shiftdim(chi2inv([1-alpha/2],dof1),-1)];

bestfit = zeros(Nchannels,Nsamples);
for ch = 1:Nchannels,
    if any(sig(ch,:)),
        lines = find(sig(ch,:));
        for L=1:length(lines),
            bestfit(ch,1:Nsamples) = bestfit(ch,1:Nsamples)+...
                sig(ch,lines(L)).*2.*real(repmat(Cstat(ch,L),[1,Nsamples]).*exp(-2*pi*complex(0,1).*repmat(harmonics(L)',[1,Nsamples]).*(0:Nsamples-1)./sampleFreq));
        end;
        data(ch,:) = data(ch,:) - bestfit(ch,:);
    end;
end;
% Recompute with lines removed
[jk,sigPow] = taperedSpectralEstimate(data,taper,pad,sample_step,keep,0); 
jk = jk.*sqrt(sample_step);

if Ktapers>1
    waitbar(0,w2,'Adaptive weighting of estimates for noise spectrum');
    for t = 1:size(jk,1)
        waitbar(t/(size(jk,1)),w2);
        Sk = squeeze(abs(jk(t,:,:)).^2);
        S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
        
        tol=.0005*sigPow(t)/pad; % Set tolerance to match PMTM
        mser=sigPow(t)*(1-concentration);
        
        temp=zeros(nKeep,1);
        S1=zeros(nKeep,1);
        b=(S0*ones(1,Ktapers))./(S0*concentration'+ones(nKeep,1)*mser'); % Percival & Walden p368a.
        wk=(b.^2).*(ones(nKeep,1)*concentration'); 
        while sum(abs(S0-S1)/nKeep)>tol
            b=(S0*ones(1,Ktapers))./(S0*concentration'+ones(nKeep,1)*mser'); % Percival & Walden p368a.
            wk=(b.^2).*(ones(nKeep,1)*concentration'); 
            S1=sum(wk'.*Sk')./ sum(wk'); % Percival & Walden p370.
            temp=S1'; 
            S1=S0; 
            S0=temp; 
        end
        dof2(t,:) = (2*(sum(wk,2)).^2./sum(wk.^2,2))';
        noisepsd(t,:) = 2*S0';
    end
    close(w2);
else
    noisepsd = abs(jk).^2;
    dof2 = 2;
end
q2 = [shiftdim(chi2inv([alpha/2],dof2),-1);shiftdim(chi2inv([1-alpha/2],dof2),-1)];

% subset = unique(dof2);
% for i=1:length(subset)
%     insubset = find(dof2==subset(i));
%     temp = chi2inv([alpha/2 1-alpha/2],subset(i));
%     q2(insubset,1:2) = repmat(temp,[length(insubset) 1]);
% end


out.freq_grid = freq_grid;
out.Ktapers = Ktapers;
out.bandwidth = bandwidth;
out.mixedSpectrum.psd = mixedpsd;
out.mixedSpectrum.dof = dof1;
if Nchannels == 1
    out.mixedSpectrum.chi2CI = [dof1.*mean(mixedpsd,1)./shiftdim(q1(2,:,:),1); dof1.*mean(mixedpsd,1)./shiftdim(q1(1,:,:),1)];
end
out.lineComponents.harmonics = harmonics;
out.lineComponents.Chat = Cstat;
out.lineComponents.Fstat = Fstat;
out.lineComponents.Pval = Pval;
out.lineComponents.sig = sig;
out.lineComponents.bestFit = bestfit;
out.backgroundSpectrum.residuals = data;
out.backgroundSpectrum.psd = noisepsd;
out.backgroundSpectrum.jk = jk;
out.backgroundSpectrum.dof = dof2;
if Nchannels == 1
    out.backgroundSpectrum.chi2CI = [dof2.*mean(noisepsd,1)./shiftdim(q2(2,:,:),1); dof2.*mean(noisepsd,1)./shiftdim(q2(1,:,:),1)];
end