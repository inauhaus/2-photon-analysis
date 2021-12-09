function T = TensorZscore_space(T)

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

Tvec = Tvec'; %this is the only difference from the "time version"

dimV = size(Tvec);
mu = mean(Tvec);
sig = std(Tvec);
Tvec = Tvec - (ones(dimV(1),1)*mu); %divide each time
Tvec = Tvec./(ones(dimV(1),1)*sig); %divide each time

Tvec = Tvec';


if ~vectorizedBit %if frames didn't start as a vector, put it back:
    T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
else
    T = Tvec;
end