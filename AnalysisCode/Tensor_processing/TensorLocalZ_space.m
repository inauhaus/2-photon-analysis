function TZ = TensorLocalZ_space(To,sig)


dim = size(To);

hh = fspecial('gaussian', [dim(1) dim(2)], sig);

Tmu = To; %preallocate
T2mu = To;
for i=1:dim(3)
    
    dum = To(:,:,i);
    
    T2mu(:,:,i) = ifft2(fft2(dum.^2).*abs(fft2(hh)));
    
    Tmu(:,:,i) = ifft2(fft2(dum).*abs(fft2(hh)));

end


Tsig = sqrt(T2mu - Tmu.^2);

TZ = (To-Tmu)./Tsig;

