function skewProj = getSkewProj(CH,hannW)

global ACQinfo

Fperiod = ACQinfo.msPerLine*ACQinfo.linesPerFrame; %frame period in ms

CHsfilt = tempsmooth(hann(round(hannW/Fperiod)),CH); %first smooth along dimension of taking max
% CHsfilt = zscore(CHsfilt,[],3);
% skewProj = mean(CHsfilt,3)-median(CHsfilt,3);
skewProj = skewness(CHsfilt,[],3);

% h = fspecial('gaussian',size(im),1);
% im = ifft2(fft2(im).*abs(fft2(h)));

function CH = tempsmooth(tkern,CH)

Idim = size(CH);
tkern = tkern/sum(tkern);
kern = zeros(Idim(1),Idim(2),length(tkern));
for i = 1:length(tkern)
    kern(:,:,i) = ones(Idim(1),Idim(2))*tkern(i);
end

smoother = zeros(size(CH));
smoother(:,:,1:length(tkern)) = kern;

smoother = abs(fft(smoother,[],3));
CH = ifft(fft(CH,[],3).*smoother,[],3);