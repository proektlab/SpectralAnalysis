function [tfse,v,w2]= hrTFspecAnalog(TS1,sampleFreq,Ktapers,freqsOfInterest,v,w2)
%function [tfse,v,w2]= hrTFspecAnalog(TS1,sampleFreq,Ktapers,freqsOfInterest,v,w2)
%
% TS1 = an analog matrix containing the timeseries of the data, as [mtrials nsample_points]
% sampleFreq = data sampling frequency (default = 1 KHz);
% Ktapers = # of tapers to use (default = 5);
%			Note that Ktapers determines the frequency resolution of the eventual output, because
% 				2*(N*sample_step*W)-1=Ktapers defines the maximum time-frequency limit you can get
%			As an example, with 2 seconds of data sampled at 1 KHz, using 5 tapers will yield 
%           a 1.5 Hz bandwidth
%
% v and w2 are optional arguments to save time when performing repeated 
% 		calculations using the same trial_length & sample_frequency
if nargin<2 sampleFreq = []; end;
if nargin<3 Ktapers = []; end;
if nargin<4 freqsOfInterest = []; end;
if nargin<5 v = []; end;
if nargin<6 w2 = [];end;

[nTrials,nSamples] = size(TS1);

if isempty(sampleFreq) sampleFreq = 1000; end;
if isempty(Ktapers)  Ktapers = 5; end
if isempty(freqsOfInterest) freqsOfInterest = [0 sampleFreq/2];
elseif length(freqsOfInterest)==1 freqsOfInterest = [0 freqsOfInterest]; end;

sample_step	= 1/sampleFreq; 
trials		= size(TS1,1); % Number of trials in the dataset
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

xk1=zeros(1,length(pad_fkeep2),Ktapers);


if (isempty(w2))
   disp(' Crossing Tapers');
   w2=zeros(Ktapers,Ktapers,nt1); % The cross-taper matrix
   [x,y]=meshgrid(1:Ktapers);
   w2=(reshape((v(ind,x(:)).*v(ind,y(:)))',[Ktapers,Ktapers,length(ind)])); 
   clear x y;
end
   
   
disp('  ...Removing linear trends...'); % For each channel, correct for any DC offset & linear trends
mrTS1=detrend(TS1')';

num_samples=trials;
tfse1a=zeros(num_samples,nt1,length(pad_fkeep2));

disp('  ...Begining Time-Frequency Power Calculation ... ');
xk1=taperedSpectralEstimate(mrTS1,v,pad,sample_step,pad_fkeep2,1);
disp(' Crossing kernel estimates trial-by-trial (using high-resolution method)');
%  Trial-by-trial, compute spectrograms by the Thomson high-resolution method
disp(['  for ' int2str(num_samples) ' trials.']);
H = waitbar(0,'Computing spectrogram');
for sample=1:num_samples       
   waitbar(sample/num_samples,H);
   sum_pk1=zeros(nt1,length(pad_fkeep2)); 
   for taper1=1:Ktapers
      for taper2=1:Ktapers
         bigw2=repmat(squeeze(w2(taper1,taper2,:)),[1 length(pad_fkeep2)]);
         bigx1_of_f=repmat((2*real(xk1(sample,:,taper1).*conj(xk1(sample,:,taper2)))),[nt1 1]);
         sum_pk1=sum_pk1 + (bigw2.*bigx1_of_f);
      end
   end
   tfse1a(sample,:,:)=sum_pk1./Ktapers;
end
clear bigx1_of_f; 
close(H);

tfse.xk1 = xk1;
tfse.bandwidth = bandwidth;
tfse.tfse=tfse1a;
tfse.freq_grid=use_freqs;
tfse.time_grid=ind*sample_step;