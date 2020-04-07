function [spect,pow] = taperedSpectralEstimate(data,taper,pad,sample_step,keep,normalize)
% function [spect,pow] = taperedSpectralEstimate(data,taper,pad,sample_step,keep,normalize)
%
% This function computes a tapered estimate of the frequency content in
% time series data.  This function implements a vectorized algorithm to
% compute the estimate, but reverts to a loop method when the size of the
% dataset is large (note: this should probably be refined).  The taper is
% passed in as a matrix, and premultiplies the data before the fft.  This
% allows apodization of the fft, to optimize for spectral leakage, etc.  A
% good choice for the taper would be the slepian sequences, computed by the
% dpss function.  Note that you can have multiple estimates computed
% simultaneously with different tapers.
%
% INPUTS:
%   data: matrix of [trials samples]
%   taper: matrix of [samples tapers] (note that this is the default output
%       format from the dpss function)
%   pad: the length to pad the fft; defaults to the next power of 2 greater
%       than the number of sample points
%   sample_step: the time in seconds between sample steps; default is 0.001
%   keep: which frequencies to keep from the pad frequencies; the default 
%       includes all frequencies up to the Nyquist (i.e., up to pad/2)
%   normalize: boolean to indicate whether the spectral estimate should be
%       normalized by sqrt(Nyquist), to give a spectrum computed from
%       spect.*conj(spect) units of power spectral density. Default is to
%       normalize, if a sample step is specified.
%
% OUTPUTS:
%  spect: a matrix containing the spectral estimates, arranged as
%       [trial,frequency,taper]
%  power: a vector containing the total power in each trial
%
% Created by Andrew Hudson, 10/12/2003.  Last updated 11/10/2004.

memlimit = 8*8*16000000; % Pick a reasonable memory limit to minimize swapping; 
% if your computer is using the scratch disk a lot
% during the execution of this function, decrease
% memlimit
if nargin<3, pad =[]; end;
if nargin<4, sample_step = []; end;
if nargin<5, keep = []; end;
if nargin<6, normalize = []; end;

if isempty(pad), pad = max(min(2^nextpow2(size(data,2)),2^17),size(data,2)); end
if isempty(sample_step), sample_step=0.5; end; % don't normalize for spectral density without a sample frequency
if isempty(keep), keep = 1:pad/2; end
if isempty(normalize), normalize = 1; end

[Ntrials,Nsamples]=size(data);
[Ksamples,Ktapers]=size(taper);

if Ksamples~=Nsamples
    if Ntrials==Ksamples
        warning(' The number of samples for data and tapers do not match; taking transpose of the data matrix to make them correspond');
        data = data';
        [Ntrials,Nsamples]=size(data);
    else
        error(' The taper must have as many points as the data');
    end
end
spect=zeros(Ntrials,length(keep),Ktapers);

if ((Ntrials*pad*Ktapers)<(memlimit/8))
    temp = repmat(data,[1 1 Ktapers]).*repmat(shiftdim(taper,-1),[Ntrials 1 1]);
    bigSpect = fft(temp,pad,2);
    spect = bigSpect(:,keep,:);
    clear temp;
elseif (Ktapers<Ntrials)&&((Ntrials*pad)<(memlimit/8)) % step down the speed for the sake of memory
    H = waitbar(0,' Computing spectral estimate');
    for t=1:Ktapers
        waitbar(t/Ktapers,H);
        temp = data.*repmat(shiftdim(taper(:,t),-1),[Ntrials 1]);
        xft1=fft(temp,pad,2);
        spect(:,:,t)=xft1(:,keep);
        clear temp;
    end
    close(H);
else
    H = waitbar(0,' Computing spectral estimate');
    for trial=1:Ntrials % Create spectral estimates, trial-by-trial
        if (mod(trial,10)==0) 
            waitbar(trial/Ntrials,H); 
        end;
        X1=data(trial,:)';
        temp = repmat(X1(1:Nsamples),1,Ktapers).*taper;
        xft1=fft(temp,pad);	
        spect(trial,:,:)=xft1(keep,:);
        clear temp;
    end
    close(H);
end

if normalize
    spect = spect.*sqrt(2*sample_step); % units of spectral density
end

if nargout>1 
    pow = var(data,1,2);
end