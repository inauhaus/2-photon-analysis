function T = TensorFilter_time(T,Lsig,Hsig)


if Lsig>0 || ~isinf(Hsig)
    
    vectorizedBit = 1;
    dim = size(T);
    
    if length(dim) == 3
        vectorizedBit = 0;
        Tvec = reshape(T,[dim(1)*dim(2) dim(3)])'; %vectorize
    else
        Tvec = T;
    end
    
    dimV = size(Tvec);
    
    if Lsig>0 && ~isinf(Hsig)
        hL = fspecial('gaussian',[dimV(1) 1],Lsig);
        hH = fspecial('gaussian',[dimV(1) 1],Hsig);
        h = hL-hH;
    elseif Lsig>0 && isinf(Hsig)
        h = fspecial('gaussian',[dimV(1) 1],Lsig);
    else 
        hL = fspecial('gaussian',[dimV(1) 1],.2); %hack to create an imppulse
        hL = hL/sum(hL);
        hH = fspecial('gaussian',[dimV(1) 1],Hsig);
        h = hL-hH;
    end

    h = h*ones(1,dimV(2));
    Tvec = ifft(abs(fft(h)).*fft(Tvec));
    
    if ~vectorizedBit %if frames didn't start as a vector, put it back to 3D:
        T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
    else
        T = Tvec;
    end

end