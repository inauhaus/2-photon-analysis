function [idy idx RFshift] = getGaborCenter(RF,sig)     

magim = RF.^2;
hh = fspecial('gaussian',size(magim),sig);
LocRMS = sqrt(ifft2(fft2(magim).*abs(fft2(hh))));

[idy idx] = find(LocRMS == max(LocRMS(:)));
dim = size(LocRMS);

LocRMS = circshift(LocRMS,[dim(1)/2-idy dim(2)/2-idx]);
RFshift = circshift(RF,[dim(1)/2-idy dim(2)/2-idx]);
RFshift = RFshift.*(hann(dim(1))*hann(dim(2))');