function aligned = removeLineNoise(aligned,sampleFreq)
%function aligned = removeLineNoise(aligned,sampleFreq);

nyquist = sampleFreq/2;

if iscell(aligned)
    matform = cell2mat(aligned)';
else
    matform = aligned;
end
[Ntrials,Nsamples] = size(matform);

oldResid = matform;

% for i = 1:floor(nyquist/60)
%     disp(['Removing any lines at ' int2str(60*i) ' Hz...']);
%     [newResid,f0] = spectralLineFilter(oldResid,sampleFreq,11,[-1 1]+60*i,2.^(nextpow2(Nsamples)+5),[],0.025); 
%     oldResid = newResid;
% end
newResid = spectralCombFilter(oldResid,sampleFreq,5,60:60:nyquist,1,2.^(min(nextpow2(Nsamples)+4,17)),[],.1);%0.025);


if all(isinf(oldResid(:)))
    oldResid = matform;
end

if iscell(aligned),
    blank = zeros(length(aligned),1);
    for i = 1:length(aligned)
        blank(i) = isempty(aligned{i});
    end
    replace = find(~blank);
    for i = 1:size(oldResid,1)
        aligned{replace(i)} = newResid(i,:)';
    end
else
    aligned = newResid;
end