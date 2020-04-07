function [Fhat,Chat,freq_grid,dof,Fs,Cs,taper] = MSMTfspectrum(data,sample_step,segmentLength,W,pad,freqsOfInterest,taper)
% function [Fhat,Chat,freq_grid,dof,Fs,Cs,taper] = MSMTfspectrum(data,sample_step,segmentLength,W,pad,freqsOfInterest,taper)
%
% 
%
%       Fhat, the fstatistics of the spectrum.  Asymptotically
%           distributed as F_2_,_{dof-2}.
%       cs, the complex spectrum; can be used to reconstruct the line
%           elements
%       freq_grid, a vector of the frequencies contained in the spectral
%           quantities

[Nchannels,Nsamples]=size(data);
if nargin<2, sample_step = []; end;
if nargin<3, segmentLength = []; end; 
if nargin<4, W = []; end;
if nargin<5, pad = []; end;
if nargin<6, freqsOfInterest = []; end;
if nargin<7, taper = [];end;

if isempty(sample_step), sample_step = 0.001; end; 
if isempty(segmentLength), segmentLength = Nsamples*sample_step; end;
Bsamples = round(segmentLength/sample_step);
Bsegments = max(ceil(2*(Nsamples-Bsamples)./Bsamples)+1,1);
Bstep = (Nsamples-Bsamples)/max(Bsegments-1,1);
Boffsets = round(Bstep*(0:Bsegments-1));
Toffsets = Boffsets*sample_step; % offset in time for phase correction

if isempty(W), W = 4/segmentLength; end;
NW = segmentLength*W;
Ktapers = ceil(2*NW) - 1;
if isempty(taper), [taper]=dpss(Bsamples,NW); taper=taper(:,1:Ktapers); end;

Nestimates = Bsegments*Ktapers;
K = zeros(Nsamples,Nestimates);
for b = 1:Bsegments,
    K(Boffsets(b)+(1:Bsamples),(b-1)*Ktapers+(1:Ktapers)) = taper;
end
Svar = norm(K'*K,'fro').^2;
clear k;
dof = 2*(Nestimates^2)/Svar;


if isempty(pad),  pad = max(2.^min(max((nextpow2(Bsamples)+4),12),16),Nsamples); end
if isempty(freqsOfInterest), freqsOfInterest = [0 1/(2*sample_step)]; 
elseif (length(freqsOfInterest)==1), freqsOfInterest = [0 freqsOfInterest]; end
freqsOfInterest = sort(freqsOfInterest);
pad_freqs=(0:pad-1)/(pad*sample_step);
keep = find(pad_freqs>=freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));
freq_grid = pad_freqs(keep);
Nkeep = length(keep);

data = detrend(data')';
segmentedData = zeros(Nchannels*Bsegments,Bsamples);
d = zeros(size(data,1),max(Bsamples,size(data,2))); 
d(1:length(data)) = data;
data = d; clear d;
phaseCorrection = zeros(Nchannels*Bsegments,length(freq_grid));
for b=1:Bsegments
    segmentedData((b-1)*Nchannels+(1:Nchannels),1:Bsamples) = data(1:Nchannels,(1:Bsamples)+Boffsets(b));
    phaseCorrection((b-1)*Nchannels+(1:Nchannels),:) = repmat(exp(-2*pi*complex(0,1)*Toffsets(b).*freq_grid),[Nchannels 1]);
end
segmentedData = detrend(segmentedData')'; % Detrend after windowing the data

% The large frequency grid sometimes causes memory problems.
brkpts = 1:1000:size(segmentedData,1);
if ~(brkpts(end) == size(segmentedData,1))
    brkpts = [brkpts size(segmentedData,1)];
elseif (length(brkpts)==1)
    brkpts(2) = 1;
end

xk(Bsegments,Nkeep,Ktapers) = 0;
for i=1:length(brkpts)-1,
    xk(brkpts(i):brkpts(i+1),:,1:Ktapers) = taperedSpectralEstimate(segmentedData(brkpts(i):brkpts(i+1),:),taper,pad,sample_step,keep,0);
end;
xk = xk.*repmat(phaseCorrection,[1 1 Ktapers]);
clear phaseCorrection;
taperPower = sum(taper(:,1:2:Ktapers),1); % even tapers are nearly balanced, so we only worry about odd ones
totalTaperPower = sum(taperPower.^2,2);
bandPower = squeeze(sum(abs(xk.^2),3))/totalTaperPower;

for ch = 1:Nchannels*Bsegments
    if (Ktapers>2)
        cs(ch,:) = (squeeze(xk(ch,:,1:2:Ktapers))*taperPower'/totalTaperPower)';
    else
        cs(ch,:) = xk(ch,:,1)*taperPower'/totalTaperPower;
    end
end

linePower = abs(cs).^2;
fs = (Ktapers-1)*linePower./max(bandPower-linePower,eps);

Cs = zeros(Nchannels,Bsegments,Nkeep);
Fs = Cs;
for b=1:Bsegments
    Cs(:,b,:) = cs((b-1)*Nchannels+(1:Nchannels),:);
    Fs(:,b,:) = fs((b-1)*Nchannels+(1:Nchannels),:);
end

clear cs fs;
Chat = zeros(Nchannels,Nkeep);
Fhat = Chat;
for n = 1:Nchannels,
    if Bsegments  > 4
        Chat(n,:) = robustLocation(squeeze(Cs(n,:,:)),'rayl',1);
        Fhat(n,:) = robustLocation(squeeze(Fs(n,:,:)),'f',2,2*Ktapers-2);
    elseif Bsegments > 1
        Chat(n,:) = mean(squeeze(Cs(n,:,:)));
        Fhat(n,:) = mean(squeeze(Fs(n,:,:)));
    else
        Chat(n,:) = (squeeze(Cs(n,1,:)));
        Fhat(n,:) = (squeeze(Fs(n,1,:)));
    end
end
end