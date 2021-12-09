function im = localXCorr2(T,sig)

%2 does linear spatial filtering of the log.

%input is tensor(x,y,T), output is image(x,y)

global ACQinfo

dim = size(T);
Tvec = reshape(T,[dim(1)*dim(2) dim(3)])';
dimV = size(Tvec);
Tvec = Tvec./(ones(dimV(1),1)*sqrt(sum(Tvec.*Tvec)));

T = reshape(Tvec',[dim(1) dim(2) dim(3)]);

%First smooth
hh=fspecial('gaussian',[dim(1) dim(2)],1);
hh = abs(fft2(hh));
for i = 1:dim(3)    
    T(:,:,i) = ifft2(fft2(T(:,:,i)).*hh);    
end

T = T-min(T(:))+.5;
T = log(T);

hh=fspecial('gaussian',[dim(1) dim(2)],sig);
hh=hh/max(hh(:));
hh = (sign(hh-.5)+1)/2;
hh = abs(fft2(hh));
for i = 1:dim(3)    
    T(:,:,i) = ifft2(fft2(T(:,:,i)).*hh);    
end
T = exp(T);
im = sum(T,3);

im = im/std(im(:));
im = im-median(im(:)); 
figure,imagesc(im)
