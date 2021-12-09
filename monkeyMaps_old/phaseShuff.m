function imShuff = phaseShuff(im) 

imMag = abs(fft2(im));
imAng = angle(fft2(im));

[M N] = size(im);

randphase = angle(fft2(rand(M,N)));
imAng = imAng+randphase;

imShuff = ifft2(imMag.*exp(1i*imAng));
