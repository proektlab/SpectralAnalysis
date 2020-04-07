function [xk,dof,biasIndicator] = suppressEigencoefficientSidelobes(yk,power,taper,concentration,pad);
% function [xk,dof,biasIndicator] = suppressEigencoefficientSidelobes(yk,power,taper,concentration,pad);
%
% This function performs the iterative sidelobe suppression outlined in
% Thomson's seminal 1982 paper.  This involves the non-linear least-squares
% solution of an integral equation to suppress the expected out-of-band 
% bias in higher order tapers.  Note that this process is data driven, and 
% hence adaptive.
%
% See either Thomson (1982) or Percival and Walden (1993) for more details.
%
% INPUTS:
%       yk, an [NtimeSegments, NfrequencyPoints] set of Fourier transforms
%              tapered by Slepian functions.
%       power, a NtimeSegments element vector specifying the total power in
%              the corresponding data segment, i.e., var(xt).
%       taper, the [NfrequencyPoints, Ktapers] matrix of the Slepian
%              functions used to create yk.
%       concentration, [Ktapers 1] vector of eigenvalues for the tapers,
%              indicating the fraction of power found within the central
%              bandwidth of the corresponding taper.  The entries of yk,
%              taper, and concentration are assumed to be ordered in
%              decreasing concentration (the default).
%
% OUTPUTS:
%       xk, the sidelobe-suppressed eigencoefficients.  Averaging their
%              magnitude is the adaptive power spectral density estimate:
%                       PSD = 2*mean(abs(xk).^2,3)/sampleFrequency
%       dof, the approximate number of degrees of freedom for each estimate
%       biasIndicator, the average over frequencies of dof/(2*Ktapers). If
%              significantly less than 1, either increase the smoothing
%              bandwidth, or apply additional prewhitening.
%
% (c) 7/27/2006, Andrew Hudson
if nargin<5, pad = []; end;
[Bsegments,Nsamples,Ktapers] = size(yk);
if isempty(pad), pad = Nsamples; end;
if (Ktapers>1) & ~isempty(concentration)
    w2 = waitbar(0,'Adaptive weighting of estimates');
    concDims = size(concentration);
    if concDims(1)>concDims(2), concentration = concentration'; end;
    if length(concentration)>Ktapers, concentration = concentration(1:Ktapers); end;
    for t = 1:Bsegments
        waitbar(t/Bsegments,w2);
        Sk = squeeze(abs(yk(t,:,:)).^2);
        
        S0=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate
        
        tol=.0005*power(t)/pad; % Set tolerance to match PMTM
        
        expectedBias = power(t)*(1-concentration); % expected out-of-band power
        
        temp=zeros(Nsamples,1);
        S1=zeros(Nsamples,1);
        dk = ones(Nsamples,Ktapers);
        while sum(abs(S0-S1)/Nsamples)>tol
            dk  = repmat(sqrt(concentration),[Nsamples 1]).*repmat(S0,[1 Ktapers])... % Thomson (1980) Eq 5.2
                ./(repmat(concentration,[Nsamples 1]).*repmat(S0,[1 Ktapers]) + repmat(expectedBias,[Nsamples 1]));
            S1 = S0;
            S0=(sum(abs(dk).^2.*Sk,2)./ sum(abs(dk).^2,2)); % Thomson (1980) Eq 5.3
        end
        xk(t,:,:) = dk.*squeeze(yk(t,:,:));
        dof(t,:) = 2*sum(abs(dk).^2,2);
    end
    biasIndicator = mean(dof./2/Ktapers,2);
    close(w2);
else
    disp('Need at least two tapers and their eigenvalues to do adaptive fitting; returning raw estimates.');
    xk = yk;
    dof = repmat(2*Ktapers,[Bsegments 1]);
    biasIndicator = [];
end
return;