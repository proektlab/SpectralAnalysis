function  [fmaxstat,width,at] = fspectrumHarmonicPeakDetection(fs,freq,f1,nHarmonics)
% function  [fmaxstat,width,at] = fspectrumHarmonicPeakDetection(fs,freq,f1,nHarmonics);
%
% 
%

if nargin<4, nHarmonics = []; end;
if isempty(nHarmonics), nHarmonics = 10; end;

ismax = imregionalmax(fs);
ismin = imregionalmin(fs);
maxlist = find(ismax);
minlist = find(ismin);

fmaxstat = zeros(1,nHarmonics);
width = zeros(1,nHarmonics);
at = zeros(1,nHarmonics);

for h = 1:nHarmonics
    [temp,mxidx] = min(abs(freq(maxlist) - f1*h));
    [temp,mnidx] = min(abs(freq(minlist) - f1*h));
    fmaxstat(h) = fs(maxlist(mxidx));
    at(h) = freq(maxlist(mxidx));
    if freq(maxlist(mxidx))>freq(minlist(mnidx)),
        check = (f1*h>=freq(minlist(mnidx))& f1*h<=freq(minlist(mnidx+1)));
        if check
            width(h) = freq(minlist(mnidx+1))-freq(minlist(mnidx));
        else
            width(h) = freq(minlist(mnidx))-freq(minlist(mnidx-1));
            fmaxstat(h) = fs(maxlist(mxidx-1));
            at(h) = freq(mxidx-1);
        end
    else
        check = (f1*h>=freq(minlist(mnidx-1))& f1*h<=freq(minlist(mnidx)));
        if check
            width(h) = freq(minlist(mnidx))-freq(minlist(mnidx-1));
        else
            width(h) = freq(minlist(mnidx+1))-freq(minlist(mnidx));
            fmaxstat(h) = fs(maxlist(mxidx+1));
            at(h) = freq(mxidx+1);
        end
    end
end