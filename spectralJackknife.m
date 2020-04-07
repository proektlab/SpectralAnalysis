function [dbMuHat,dbSigmaHat] = spectralJackknife(data,jackOverTapers)
% function [dbMuHat,dbSigmaHat] = spectralJackknife(data,jackOverTapers)
% 
% This function computes the mean and standard deviation in db of a dataset.
% Each row of the data is assumed to be an independent trial.  The data may
% be a vector or multidimensional.  If the logical flag, jackOverTapers, is set,
% the jackknife will act along the first and last dimensions of the data.

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
if length(dims)>2
   jk = reshape((1-eye(N))*reshape(data, [N prod(dims(2:end))])/(N-1),dims);
else
   jk = (1-eye(N))*data./(N-1);
end
% Make relevant calculations
dbMuHat=10*log10(abs(mean(data,1)));
dbSigmaHat=sqrt(((N-1)/N).*sum((10*log10(abs(jk)) - repmat(dbMuHat,[N 1 1])).^2,1));

