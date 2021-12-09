function T = TensorL2Normalize_time(T)

%If T is 3 dimensions, assumes time is the 3rd
%If T is 2 dimensions (vectorized frames), assumes time is 2nd

vectorizedBit = 1;
dim = size(T);

if length(dim) == 3
    vectorizedBit = 0;
    Tvec = reshape(T,[dim(1)*dim(2) dim(3)])'; %vectorize
else
    Tvec = T;
end


dimV = size(Tvec);
L2 = sqrt(sum(Tvec.*Tvec));
Tvec = Tvec./(ones(dimV(1),1)*L2); %divide each time course by norm


if ~vectorizedBit %if frames didn't start as a vector, put it back:
    T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
else 
    T = Tvec;
end