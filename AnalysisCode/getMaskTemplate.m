function [maim sigim muim] = getMaskTemplate

%%
hh = makeMapFilter;

global f0m

f0 = zeros(size(f0m{1},1),size(f0m{1},2),length(f0m)-1);
for i = 1:length(f0m)-1
   
    f0(:,:,i) = ifft2(abs(fft2(hh)).*fft2(f0m{i}));
    
end
%%
maim = max(f0,[],3);
sigim = std(f0,[],3);
muim = mean(f0,3);

mi = prctile(maim(:),.1);
ma = prctile(maim(:),99.5);
maim(find(maim>ma)) = ma;
maim(find(maim<mi)) = mi;

mi = prctile(sigim(:),.1);
ma = prctile(sigim(:),99.5);
sigim(find(sigim>ma)) = ma;
sigim(find(sigim<mi)) = mi;

mi = prctile(muim(:),.1);
ma = prctile(muim(:),99.5);
muim(find(muim>ma)) = ma;
muim(find(muim<mi)) = mi;

figure, subplot(1,3,1)
imagesc(maim)
title('max')

subplot(1,3,2)
imagesc(sigim)
title('std dev')

subplot(1,3,3)
imagesc(muim)
title('mean')

%%

maskTemps.muim = muim;
maskTemps.sigim = sigim;
maskTemps.maim = maim;
%%

global AUE

path = 'c:\';
filename = strcat(path,AUE,'_maskTemp');
uisave('maskTemps',filename)


