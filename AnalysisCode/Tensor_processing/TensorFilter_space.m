function To = TensorFilter_space(To,LPsig,HPsig)


if LPsig>0 && ~isinf(HPsig)
    
    dim = size(To);
    if LPsig>0
        hLP = fspecial('gaussian', [dim(1) dim(2)], LPsig);
    else
        hLP = zeros(dim(1),dim(1));
        hLP(round(dim(1)/2),round(dim(2)/2)) = 1;
    end
    
    if ~isinf(HPsig)
        hHP = fspecial('gaussian', [dim(1) dim(2)], HPsig);
    else
        hHP = 0;
    end
    
    hh = hLP-hHP;
    
    
    for i=1:dim(3)
        
        dum = To(:,:,i);
        dum = ifft2(fft2(dum).*abs(fft2(hh)));
        To(:,:,i) = dum;
        
    end
    
end

