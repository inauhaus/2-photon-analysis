function rTemp = findBadframes(T)

%The goal here is to ID frames that are very different from the anatomy
%template

global maskS

%2x2 binning in space
Atemp = maskS.anatomyTemplate;
A = Atemp(1:2:end,1:2:end,:);
B = Atemp(1:2:end,2:2:end,:);
C = Atemp(2:2:end,1:2:end,:);
D = Atemp(2:2:end,2:2:end,:);
Atemp = A+B+C+D;
Atemp = Atemp(:);
Atemp = (Atemp-mean(Atemp(:)))/std(Atemp(:));

A = T(1:2:end,1:2:end,:);
B = T(1:2:end,2:2:end,:);
C = T(2:2:end,1:2:end,:);
D = T(2:2:end,2:2:end,:);
To = A+B+C+D;
dim = size(To);
To = reshape(To,[dim(1)*dim(2) dim(3)])';
dimV = size(To);


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
%ToI = To(1:2:end,:);

muTo = mean(To,2)*ones(1,dimV(2));
sigTo = std(To,[],2)*ones(1,dimV(2));
To = (To-muTo)./(sigTo);

rTemp = (To*Atemp); %Correlation with the template, at each frame

