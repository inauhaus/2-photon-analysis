function [T rTemp] = TensorNaNBadFrames(T,Temp,Zthresh,rFloor)


vectorizedBit = 1;
dim = size(T);
if length(dim) == 3
    Tvec = reshape(T,[dim(1)*dim(2) dim(3)])';   
    vectorizedBit = 0;
else
    Tvec = T;
end

Temp = Temp(:);

TvecZ = TensorZscore_space(Tvec);
Temp = (Temp-mean(Temp(:)))/std(Temp(:));

rTemp = (TvecZ*Temp)/size(TvecZ,2); %Correlation coef with the template, at each frame

%figure, plot(rTemp)

thresh = prctile(rTemp,50) - Zthresh*std(rTemp);

%thresh = 0.6;
idBadTimes = find(rTemp<thresh | rTemp<rFloor);

Tvec(idBadTimes,:) = NaN;


if ~vectorizedBit %if frames didn't start as a vector, put it back to 3D:
    T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
else 
    T = Tvec;
end