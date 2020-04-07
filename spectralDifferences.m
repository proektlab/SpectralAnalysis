function [whereSigDiff,dbDiff]=spectralDifferences(dbMuHat1,dbSigmaHat1,dbMuHat2,dbSigmaHat2,alpha);
%function [whereSigDiff,dbDiff]=spectralDifferences(dbMuHat1,dbSigmaHat1,dbMuHat2,dbSigmaHat2,alpha);
%
% This function is configured to accept vector (e.g., spectrum) or multidimensional 
%  (e.g., spectrogram) data.

if nargin<5
   p=0.975;
else
   p=1-alpha/2;
end

dbDiff=dbMuHat1-dbMuHat2; % Compute the difference in the spectra in dBs 

whereSigDiff = zeros(size(dbMuHat1));
lowerLimit1 = dbMuHat1 - dbSigmaHat1.*norminv(p,0,1);
upperLimit1 = dbMuHat1 + dbSigmaHat1.*norminv(p,0,1);
lowerLimit2 = dbMuHat2 - dbSigmaHat2.*norminv(p,0,1);
upperLimit2 = dbMuHat2 + dbSigmaHat2.*norminv(p,0,1);

whereSigDiff = (dbMuHat2<lowerLimit1 | dbMuHat2>upperLimit1) ...
   & (dbMuHat1<lowerLimit2 | dbMuHat1>upperLimit2);   

