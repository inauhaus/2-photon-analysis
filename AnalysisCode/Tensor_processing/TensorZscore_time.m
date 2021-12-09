function T = TensorZscore_time(T)

vectorizedBit = 1;
dim = size(T);
if length(dim) == 3
    Tvec = reshape(T,[dim(1)*dim(2) dim(3)])';
    vectorizedBit = 0;
else
    Tvec = T;
end

dimV = size(Tvec);

mu = mean(Tvec);
sig = std(Tvec);

Tvec = Tvec - (ones(dimV(1),1)*mu); %divide each time course by norm
Tvec = Tvec./(ones(dimV(1),1)*sig); %divide each time course by norm


if ~vectorizedBit %if frames didn't start as a vector, put it back to 3D:
    T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
else 
    T = Tvec;
end