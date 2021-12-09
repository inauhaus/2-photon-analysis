function ToI = cleanTensor(T)

%Will return an array that is downsampled 2x

A = T(1:2:end,1:2:end,:);
B = T(1:2:end,2:2:end,:);
C = T(2:2:end,1:2:end,:);
D = T(2:2:end,2:2:end,:);
To = A+B+C+D;
dim = size(To);
To = reshape(To,[dim(1)*dim(2) dim(3)])';

% h = 1./fft(mean(imR,2));
% h = h*ones(1,size(imR,2));
% imR = ifft(abs(fft(h)).*fft(imR));

%Smooth each pixel in time
tsig = 1;
if tsig>0
    h = fspecial('gaussian',[size(To,1) 1],tsig);
    h = h*ones(1,size(To,2));
    To = ifft(abs(fft(h)).*fft(To));
end
ToI = To(1:2:end,:);

% Subtract polynomial fit from each pixel (i.e. HP filter)
H = (1:size(ToI,1))';
H = [H H.^2 H.^3 ones(length(H(:,1)),1)];
slps = inv(H'*H)*H'*ToI;
imRfit = H*slps;
ToI = ToI-imRfit;


%%

%Remove first singular value

[U,S,V] = svd(ToI','econ');
C1 = U(:,1)*V(:,1)'*S(1,1);
ToI = ToI-C1';

ToI = reshape(ToI',[dim(1) dim(2) size(ToI,1)]);




% dim = size(imXdum);
% dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
% imXdum = interp1(1:dim(1),imXdum,linspace(1,dim(1),dimI(1)));
% imXdum = interp1(1:dim(2),imXdum',linspace(1,dim(2),dimI(2)))';





