function TemplatefromStimulusTriggered

%This assumes that the first dimension of kernelsIm is orientation.  It
%then sums across the other dimensions and computes the unnormalized tuning
%selectivity.

global kernelsIm

dim = size(kernelsIm);

for i = 1:dim(1)
    imT{i} = 0;
    for j = 1:dim(2)
        for k = 1:dim(3)
            for l = 1:dim(4)
               	dum = kernelsIm{i,j,k,l};
                imT{i} = imT{i}+dum;
            end
        end
    end
end

%%
Tid = 2:10;
clear im
im = zeros(size(kernelsIm{1},1),size(kernelsIm{1},2));
for i = 1:dim(1)
    im(:,:,i) = mean(imT{i}(:,:,Tid),3);
end

dori = 180/size(im,3);
oridom = 0:dori:180-dori;

angmap = 0;
for i = 1:dim(1)
    angmap = angmap + im(:,:,i)*exp(1i*2*oridom(i)*pi/180);
end

h = fspecial('gaussian',size(angmap),1);
angmaph = ifft2(abs(fft2(h)).*fft2(angmap));
orimap = angle(angmaph)/2*180/pi;
id = find(orimap<0);
orimap(id) = 180+orimap(id);

magmap = abs(angmaph);
ma = prctile(magmap(:),99.5);
mi = prctile(magmap(:),1);
magmap(find(magmap>ma)) = ma;
magmap(find(magmap<mi)) = mi;

magmap = magmap-min(magmap(:));
magmap = magmap/max(magmap(:));

figure,imagesc(orimap,'AlphaData',magmap),colormap hsv
figure,imagesc(magmap), colormap gray



%%


dim = [size(kernelsIm{1},1) size(kernelsIm{1},2)];
vim = zeros(dim(1),dim(2),2);
for i = 1:prod(size(kernelsIm))    
   kdum = kernelsIm{i};
   vim(:,:,2) = max(kdum,[],3);
   vim(:,:,1) = max(vim,[],3);    
end

im = vim(:,:,1).^2;

h = fspecial('gaussian',dim,1);
imh = ifft2(abs(fft2(h)).*fft2(im));

ma = prctile(imh(:),99.8);
mi = prctile(imh(:),.1);

imh(find(imh>ma)) = ma;
imh(find(imh<mi)) = mi;

figure,imagesc(imh), colormap gray


% maskS.im = imh