function im = localXCorr3(Traw,sig)

%2 does linear spatial filtering of the log.

%input is tensor(x,y,T), output is image(x,y)

global ACQinfo

%%
T = Traw;

idNaN = find(isnan(T(1,1,:)));
T(:,:,idNaN) = [];

%sig = 40;
%T = TensorLocalZ_space(T,sig);

mag = ACQinfo.SBInfo.config.magnification;

%T = TensorFilter_time(T,.3,10);

dim = size(T);

%T = TensorL2Normalize_time(T);

%First get local neighborhood
nsig1 = 3*mag;
hh1=fspecial('gaussian',[dim(1) dim(2)],nsig1);
nsig2 = 20*mag;
hh2=fspecial('gaussian',[dim(1) dim(2)],nsig2);
hh = hh1-hh2;

%hh = hh1;

hh(find(hh == max(hh(:)))) = 0;
hh = abs(fft2(hh));
Tn = T;
for i = 1:size(Tn,3)  
    Tn(:,:,i) = ifft2(fft2(Tn(:,:,i)).*hh);    
end

sig = 10;
%Tn = TensorLocalZ_space(Tn,sig);
%T = TensorLocalZ_space(T,sig);

Tn = TensorZscore_time(Tn);
T = TensorZscore_time(T);

%Tn = TensorL2Normalize_time(Tn);
%T = TensorL2Normalize_time(T);


%Tn = getTcourseNorm(Tn);


%im = localSVD(T,3);

im = mean(Tn.*T,3);

im = im-prctile(im(:),70);
im = phi(im);
%im = medfilt2(im,[3 3]);
% 

% im = im/std(im(:));
%im = phi(im-prctile(im(:),5)).^2; 

ma = prctile(im(:),99);
mi = prctile(im(:),1);

figure,imagesc(im,[mi ma]), colormap gray
%%

function Tnorm = getTcourseNorm(T)

dim = size(T);
Tvec = reshape(T,[dim(1)*dim(2) dim(3)])';
dimV = size(Tvec);
nom = sqrt(sum(Tvec.*Tvec));
Tvec = Tvec./(ones(dimV(1),1)*nom); %divide each time course by norm

Tnorm = reshape(Tvec',[dim(1) dim(2) dim(3)]); %reshape


