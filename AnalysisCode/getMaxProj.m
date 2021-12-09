function maxProj = getMaxProj(CH,hannW)

global ACQinfo

Fperiod = ACQinfo.msPerLine*ACQinfo.linesPerFrame; %frame period in ms

CHsfilt = tempsmooth(hann(round(hannW/Fperiod)),CH); %first smooth along dimension of taking max

maxProj = max(CHsfilt,[],3);


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