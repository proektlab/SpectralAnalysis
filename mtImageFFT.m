function [Out] = mtImageFFT(image,tapers, imDim)
% As input accepts an image (2D array). tapers is the number of tapers.
% Will assume the optimal timebadnwidth as implemented in swTFspecAnalog.
% imDIm is a two element vector that gives the size of the image in units
% of space. If that is not provided will assume dimensions of [1 1]. 
% Outputs: out.freq1=frequency along first dimension, out.freq2=frequency
% along second dimension. out.spectrum is the squared modulus of 2D FFt.
% out.phase is the spatial phase. 
if nargin==2
    imDim=[1,1];
end

[s1, s2]=size(image);               % dimensions of the image in units of pixels
[out,~,~] = swTFspecAnalog(image,s2,tapers,[],s2,1,[],[],[],[],[]); % perform multitaper along second dimension
R=squeeze(out.tfse);            % get the absolute value of the fft
Ph=squeeze(out.phase);          % get the phase
[X, Y]=pol2cart(Ph, sqrt(R));   % convert to polar coordinates;
CompSpect=X+1i.*Y;              % reconstruct the complex valued spectrum
freq2=out.freq_grid;            % extract frequency
freq2=freq2./imDim(2);          % convert to units of space if supplied

% now repeat the spectral analysis on in the other direction
[out,taper,concentration] = swTFspecAnalog(CompSpect.',s1,tapers,[],s1,1,[],[],[],[],[]);
% repeat the conversion
R=squeeze(out.tfse);
Ph=squeeze(out.phase);
[X, Y]=pol2cart(Ph, sqrt(R));
CompSpect2=X+1i.*Y;
CompSpect2=CompSpect2.';
freq1=out.freq_grid;
freq1=freq1./imDim(1);

Out.freq1=freq1;
Out.freq2=freq2;
Out.phase=angle(CompSpect2);
Out.spectrum=abs(CompSpect2).^2;


end

