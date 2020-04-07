function [artifacts,pvals,epvals] = lfpArtifactDetection(aligned,doPlot);
% function [artifacts,pvals,epvals] = lfpArtifactDetection(aligned,doPlot);
%
% lfpArtifactDetection is useful for analyzing LFP data recordings that
% have been aligned to some behavioral event in order to prevent electrical
% (and possibly mechanical) artifacts from contaminating the analysis of
% the data.
%
% Given an alignment structure (cell array containing column vectors of
% data) or a raster matrix (matrix of data arranged as [trials timepoints]),
% lfpArtifact rejection attempts to model the variance within each trial
% with a gamma distribution, and to flag those trials that have a p value
% for the fitted distribution greater than 1-1/Ntrials or 0.95, whichever
% is greater.  A similar screening is made based upon the maximum excursion
% from zero that occurs during the trial.
%
% INPUTS:
%   aligned, either a cell array where each cell contains a column vector
%           of data for a trial, or a data matrix where each row contains a
%           trial.
% OUTPUTS:
%   artifacts, a list of the trials that appear to be artifacts
%   ppvals, a vector of approximate p values for the power in each trial.
%   epvals, a vector of approximate p values for the extrema in each trial.
%
% Version 0.1 
% ***PLEASE NOTE THAT THIS IS STILL A BETA VERSION AND THE FINAL
%   VERSION OF THIS CODE MAY BE ALTERED***
% Written by Andrew Hudson, 2/7/2005. 

if nargin<2, doPlot = []; end;
if isempty(doPlot), doPlot = 0; end;

if iscell(aligned)
    matform = cell2mat(aligned)';
else
    matform = aligned;
end
N = size(matform,1);

if N>4
    pow = var(matform');
    pow(find(pow==0))=eps;
    %powcutoff = exp(1.5*iqr(log(pow)) + median(log(pow))); % make the fitting robust to outliers
    powcutoff = 1.5*iqr(pow) + median(pow); % make the fitting robust to outliers
    bpow = gamfit(pow(find(pow<powcutoff)));
    powp = gamcdf(pow,bpow(1),bpow(2));
    powarts = find(powp>max(1-1/N,0.99));
    
    maxes = max((matform.^2)');
    maxcutoff = exp(1.5*iqr(log(maxes)) + median(log(maxes))); % or 3* the linear scaling
    %maxcutoff = 3*iqr(maxes) + median(maxes); % or 3* the linear scaling
    bmax = gamfit(maxes(find(maxes<maxcutoff)));
    maxp = gamcdf(maxes,bmax(1),bmax(2));
    maxarts = find(maxp>max(1-1/N,0.99));
    
    if doPlot
        figure;
        pcolor(matform);
        colormap('gray');
        shading('flat');
        hold on;
        xlim = get(gca,'xlim');
        plot(repmat(xlim,[length(powarts) 1])',repmat(powarts'+0.5,[1 2])','r');
        plot(repmat(xlim,[length(maxarts) 1])',repmat(maxarts'+0.5,[1 2])','c');
        plot(repmat(xlim,[length(intersect(maxarts,powarts)) 1])',repmat(intersect(powarts,maxarts)'+0.5,[1 2])','m');
    end
    
    artifacts = union(powarts,maxarts);
    pvals = powp;  % Note that the power p values are returned,
                   % but both power and extrema figure into the artifact
                   % labelling.  This is because the power and the extrema
                   % are not totally independent, so you can't just
                   % multiply the probabilities to find the p value of that
                   % combination, so I'm wussing out and only returning one
                   % p value.
else
    artifacts = [];
    pvals = NaN;
    epvals = NaN;
end