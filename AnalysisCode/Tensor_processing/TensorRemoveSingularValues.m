function T = TensorRemoveSingularValues(T,N)

%Remove N singular value

vectorizedBit = 1;
dim = size(T);

if length(dim) == 3
    vectorizedBit = 0;
    Tvec = reshape(T,[dim(1)*dim(2) dim(3)])'; %vectorize
else
    Tvec = T;
end

dimV = size(Tvec);
    

medmat = ones(size(Tvec,1),1)*mean(Tvec); %Center of "cloud"

Tvec = Tvec - medmat; %Center the "cloud"
%T = T./(ones(size(T,1),1)*std(T)); %Center the "cloud"


[U,S,V] = svd(Tvec','econ');

C = 0;
for i = 1:N   
    C = C + U(:,i)*V(:,i)'*S(i,i);
end

%C = U(:,1)*V(:,1)'*S(1,1) ;
%C = C + U(:,2)*V(:,2)'*S(2,2) + U(:,3)*V(:,3)'*S(3,3) ;

Tvec = Tvec-C';

Tvec = Tvec + medmat; %Put the cloud back

if ~vectorizedBit %if frames didn't start as a vector, put it back to 3D:
    T = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape
else
    T = Tvec;
end
