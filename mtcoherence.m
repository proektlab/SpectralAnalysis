function [output,v]= mtcoherence(TS1,TS2,sample_step,Ktapers,v)
%function [output,v]= mtcoherence(TS1,TS2,sample_step,Ktapers,v)
%
% INPUTS:
%     TS1, a matrix of time series, arranged as [trials, sample_points]
%     TS2, a matrix of time series, arranged as [trials, sample_points]
%     sample_step, a scalar specifying the amount of time between sample
%           points, in seconds; default: 0.005, or 200 Hz
%     Ktapers, a scalar indicating how many slepian tapers to include
%           default: 5
%     v, an optional input, containing precomputed slepian tapers
%
% OUTPUS:
%     output, a structure with the following fields:
%     psd1 & psd2, an estimate of the psd for each of the time series
%		psd1err & psd2err, an estimate of their jackknife standard error
% 		coh, an estimate of the coherency between the signals
% 		cohCI, a jackknife estimate of 99% confidence limits for the coherence
%		whereSigUP & sigCoh, uniform phase distribution estimates of the
%		    significant coherence, at 99%
%     frequencies, a list of the frequencies for which the psd is
%         provided
%     v, as above

[trials,trial_length]=size(TS1);
if trial_length == 1
   TS1=TS1';
   TS2=TS2';
   trial_length = trials;
   trials = 1;
end
TS1 = detrend(TS1')';
TS2 = detrend(TS2')';

clear output input;
if nargin<3
   sample_step = 0.005; % default to 200 Hz
end;
if nargin<4
   Ktapers = 5;
end

bandwidth= (Ktapers+1)/(trial_length*sample_step); % The limit on frequency resolution is governed by the number of tapers
fs			=1/(trial_length*sample_step); % Fundamental frequency in Hz
total_harms	=(1/(2*sample_step))/fs; % Nyquist divided by frequency of 1st harmonic (fundamental frequency)
all_freqs=fs*(1:total_harms); % Define frequency grid
hgrid		=1:total_harms; % Points in frequency grid
freqs_on_hgrid	=fs:bandwidth:(total_harms*fs); % Sampling determined by frequency grid
NW			=(Ktapers+1)/2; % Time-frequency product defines the limit on estimates
pad		=2.^nextpow2(16*trial_length); % FFT runs fastest when a power of 2
ki			=[0 freqs_on_hgrid max(freqs_on_hgrid)+bandwidth];
ka			=ki/fs;
pad_freqs=interp1(ka,ki,fs:(total_harms/(pad/2)):total_harms);

if (total_harms*fs)<200
   pad_fkeep2 = 1:ceil(length(pad_freqs)/240):length(pad_freqs); 
else
   temp = min(find(pad_freqs>=200));
   pad_fkeep2 = 1:ceil(temp/240):temp;
   clear temp;
end;
use_freqs = pad_freqs(pad_fkeep2);

if ((nargin<5)|(isempty(v)))
   v=dpss(trial_length,NW);			%generate the tapers
end
v=v(:,1:Ktapers);

rand('seed',sum(clock));
spec1 = zeros(trials,length(pad_fkeep2),Ktapers);
spec2 = spec1;
jks11 = spec1; jks22 = spec1; jks12 = spec1;
disp('Begining Power Calculation');
spec1 = taperedSpectralEstimate(TS1,v,pad,sample_step,pad_fkeep2);
spec2 = taperedSpectralEstimate(TS2,v,pad,sample_step,pad_fkeep2);

s11 = spec1.*conj(spec1);
s22 = spec2.*conj(spec2);
s12 = spec1.*conj(spec2);
disp('Determining confidence limits for coherence');
[cohCI,phaseCI] = jackCoherenceCI(s11,s22,s12,0.01);

disp('Testing for significance by uniform phase');
dims=size(spec1);
cohMat = s12./sqrt(s11.*s22);
coh = mean(mean(s12,3),1)./sqrt(mean(mean(s11,1),3).*mean(mean(s22,1),3));
if trials > 1
   [whereSigUP,sigCohUP,uniformPhaseP]=uniformPhaseTest(mean(cohMat,3),0.01);
else
   [whereSigUP,sigCohUP,uniformPhaseP]=uniformPhaseTest(permute(squeeze(cohMat),[2 1]),0.01);
end
output.psd1=s11;
output.psd2=s22;
output.coh=coh;
output.cohCI=cohCI;
output.phaseCI = phaseCI;
output.uniformPhaseP=uniformPhaseP;
output.cohMat=cohMat;
output.whereSigUP=whereSigUP;
output.sigCoh=coh.*whereSigUP;
output.frequencies = use_freqs;

%%%%
function [magCI,phaseCI] = jackCoherenceCI(s11,s22,s12,criterion);
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

jkE = jkC./abs(jkC);   % phase factor
jkMeanE = mean(jkE,1); 
jkStdErrPhi = sqrt(2*(N-1)*(1-abs(jkMeanE))); % Thomson and Chave (1991), p91, eq 2.63
%meanPhi = angle(mean(s11,1)./sqrt(mean(s11,1).*mean(s22,1)));
%phaseCI = [meanPhi + CInum.*jkStdErrPhi; meanPhi - CInum.*jkStdErrPhi];
phaseCI = [CInum.*jkStdErrPhi; - CInum.*jkStdErrPhi];

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
