function out=circavg(inim,nstep);

% function out=circavg(inim,nstep);
% 
%CIRCAVG computes the circular average of an image INIM. Most used to compute
%the circular average of a 2D fft amplitude spectrum
%abs(fftshift(fft2(inim)). 
%NSTEP specifies the number of steps along which the average is going to be
%made. NSTEP shouldn't be higher than half the size of the smaller side of the
%image array (which is the default value if NSTEP is omitted).
%
%The output is a vector of size 1 x NSTEP.
%
%written by C. Jacques September 2010

%check inputs
if ndims(inim)>2;
    error('input should be a 2D array');
end

[H,W]=size(inim);
if nargin==1; nstep=floor(min(H,W)/2);end

if nstep>floor(min(H,W)/2);
    error('number of steps should not be higher than niquist (image-size/2)');
end


[f1,f2] = freqspace([H W],'meshgrid');

r=sqrt(f1.^2 + f2.^2);

ind=round(r*nstep)+1;

for ii=1:nstep
    indii=find(ind==ii);
    out(ii)=mean(inim(indii));
end

end