function [tfse,v,w2]= hrTFcoherence(TS1,TS2,sampleFreq,Ktapers,freqsOfInterest,v,w2)
%function [tfse,v,w2]= hrTFcoherence(TS1,sampleFreq,Ktapers,freqsOfInterest,v,w2)
%
% TS1 = an analog matrix containing the timeseries of the data, as [mtrials nsample_points]
% sample_step = time (s) between samples, i.e., inverse of data sampling frequency (default = 0.001);
% Ktapers = # of tapers to use (default = 5);
%			Note that Ktapers determines the frequency resolution of the eventual output, because
% 				2*(N*sample_step*W)-1=Ktapers defines the maximum time-frequency limit you can get
%			As an example, with 2 seconds of data sampled at 1 KHz, using 5 tapers will yield 
%           a 1.5 Hz bandwidth
%
% v and w2 are optional arguments to save time when performing repeated 
% 		calculations using the same trial_length & sample_frequency

if nargin<3 sampleFreq = []; end;
if nargin<4 Ktapers = []; end;
if nargin<5 freqsOfInterest = []; end;
if nargin<6 v = []; end;
if nargin<7 w2 = [];end;

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(Ktapers)  Ktapers = 5; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;

sample_step	= 1/sampleFreq; 
T				= size(TS1,2); % Length of a trial in sample points; determines analysis window
bandwidth	= (Ktapers+1)/(T*sample_step); % The limit on frequency resolution is governed by the number of tapers
fs				=1/(T*sample_step); % Fundamental frequency in Hz
total_harms	=(1/(2*sample_step))/fs; % Nyquist divided by frequency of 1st harmonic (fundamental frequency)
all_freqs	=fs*(1:total_harms); % Define frequency grid
hgrid			=1:total_harms; % Points in frequency grid
freqs_on_hgrid	=fs:bandwidth:(total_harms*fs); % Sampling determined by frequency grid
Tstep			=0.020/sample_step; %time window stepsize is 20 msec by default 
NW				=(T*sample_step)*(bandwidth/2); % Time-frequency product defines the limit on estimates
pad			=max(2^14,2^nextpow2(T*16)); % FFT runs fastest when a power of 2

ind=[1:Tstep:T]; % Create the time grid for the spectrograms
nt1=length(ind); % Number of time points in the grid

freqsOfInterest = sort(freqsOfInterest);
pad_freqs = sampleFreq*(1:pad)/pad;
pad_fkeep2 = find(pad_freqs>freqsOfInterest(1)&pad_freqs<=freqsOfInterest(2));
frqStep = round(length(pad_fkeep2)/120);
pad_fkeep2 = pad_fkeep2(1):frqStep:pad_fkeep2(end);

use_freqs=pad_freqs(pad_fkeep2);


if isempty(v)
   disp(' Calculating tapers');
   v=dpss(T,NW);			%generate the tapers
end
v=v(:,1:Ktapers);


ind=[1:Tstep:T]; % Create the time grid for the spectrograms
nt1=length(ind); % Number of time points in the grid
if (nargin<6)
   w2=[];
end
if (isempty(w2))
   disp(' Crossing Tapers');
   w2=zeros(Ktapers,Ktapers,nt1); % The cross-taper matrix
   [x,y]=meshgrid(1:Ktapers);
   w2=(reshape((v(ind,x(:)).*v(ind,y(:)))',[Ktapers,Ktapers,length(ind)])); 
   clear x y;
end
   
   
disp('  ...Removing DC-offset...');
% For each channel, correct for any DC offset 
mrTS1=TS1-repmat(mean(TS1,2),[1 size(TS1,2)]);
mrTS2=TS2-repmat(mean(TS2,2),[1 size(TS2,2)]);

tfse1a=zeros(nTrials,nt1,length(pad_fkeep2));
tfse2a=zeros(nTrials,nt1,length(pad_fkeep2));

xk1=zeros(1,length(pad_fkeep2),Ktapers);
xk2=xk1;

disp('  ...Begining Time-Frequency Power Calculation ... ');
xk1=taperedSpectralEstimate(mrTS1,v,pad,sample_step,pad_fkeep2);
xk2=taperedSpectralEstimate(mrTS2,v,pad,sample_step,pad_fkeep2);

disp(' Crossing kernel estimates trial-by-trial (using high-resolution method)');
%  Trial-by-trial, compute spectrograms by the Thomson high-resolution method
disp(['  for ' int2str(nTrials) ' trials.']);
H = waitbar(0,'Computing spectrogram');
for trial=1:nTrials       
   %disp([ '  trial ' int2str(trial) ' of ' int2str(nTrials) ]);
   waitbar(trial/nTrials,H);
   sum_pk1=zeros(nt1,length(pad_fkeep2)); 
   sum_pk2=zeros(nt1,length(pad_fkeep2)); 
   sum_pk12=zeros(nt1,length(pad_fkeep2)); 
   for taper1=1:Ktapers
      for taper2=1:Ktapers
         bigw2=repmat(squeeze(w2(taper1,taper2,:)),[1 length(pad_fkeep2)]);
         bigx1_of_f=repmat((2*real(xk1(trial,:,taper1).*conj(xk1(trial,:,taper2)))),[nt1 1]);
         bigx2_of_f=repmat((2*real(xk2(trial,:,taper1).*conj(xk2(trial,:,taper2)))),[nt1 1]);
         sum_pk1=sum_pk1 + (bigw2.*bigx1_of_f);
         sum_pk2=sum_pk2 + (bigw2.*bigx2_of_f);
         sum_pk12=sum_pk12 + (bigw2.*repmat(xk1(trial,:,taper1).*conj(xk2(trial,:,taper2)),[nt1 1]));
      end
   end
   tfse1a(trial,:,:)=sum_pk1./Ktapers;
   tfse2a(trial,:,:)=sum_pk2./Ktapers;
   tfcoh(trial,:,:)=sum_pk12./sqrt(abs(sum_pk1).*abs(sum_pk2));
end		%for trial
clear bigx1_of_f; 
close(H);

tfse.xk1 = xk1;
tfse.xk2 = xk2;
tfse.bandwidth = bandwidth;
tfse.tfse1=tfse1a;
tfse.tfse2=tfse2a;
tfse.tfcoh=tfcoh;
tfse.freq_grid=use_freqs;
tfse.time_grid=ind*sample_step;

[whereSigUP,sigCohUP,uniformPhaseP]=uniformPhaseTest(tfse.tfcoh,0.01);
tfse.uniformPhaseP = uniformPhaseP;


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
