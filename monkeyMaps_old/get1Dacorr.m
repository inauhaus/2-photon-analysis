function [odacorr sfacorr dom odspec_pad sfspec_pad fdom_pad] = get1Dacorr(imod,imsf)

% open('E:\SFODfigs\AutoCorr\ab8_0_17_OD_2DAC.fig')
% imod = get(gcf,'Cdata');
% open('E:\SFODfigs\AutoCorr\ab8_0_17_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ab8_0_17_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ab8mats\','imod','imsf','imori')
%
% open('E:\SFODfigs\AutoCorr\ab9_0_102_OD_2DAC.fig')
% imod = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ab9_0_102_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ab9_0_102_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ab9mats\','imod','imsf','imori')
%
% open('E:\SFODfigs\AutoCorr\ac0_0_56_OD_2DAC.fig')
% imod = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac0_0_56_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac0_0_56_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ac0mats\','imod','imsf','imori')
%
% open('E:\SFODfigs\AutoCorr\ac1_0_115_OD_2DAC.fig')
% imod = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_0_115_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_0_115_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ac1mats1\','imod','imsf','imori')
%
% open('E:\SFODfigs\AutoCorr\ac1_1_11_OD_2DAC.fig')
% imod = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_1_11_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_1_11_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ac1mats2\','imod','imsf','imori')
%
% open('E:\SFODfigs\AutoCorr\ac1_2_19_OD_2DAC.fig')
% imod = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_2_19_SF_2DAC.fig')
% imsf = get(gco,'Cdata');
% open('E:\SFODfigs\AutoCorr\ac1_2_19_ORI_2DAC.fig')
% imori = get(gco,'Cdata');
% save('E:\SFODfigs\Autocorr\ac1mats3\','imod','imsf','imori')

%%

[xmicperpix ymicperpix] = getImResolution(1);
micperpix = (xmicperpix+ymicperpix)/2;

Win = hann(size(imod,1))*hann(size(imod,2))';

id = find(isnan(imod));
imod(id) = 0;
imoddum = imod-mean(imod(:));
imoddum = imoddum.*Win;
OD_Pow = fftshift(fftshift(abs(fft2(imoddum)),1),2);
%figure,imagesc(OD_Pow)

odradon = radon(imoddum,[0:10:170]);


id = find(isnan(imsf));
imsf(id) = 0;
imsfdum = imsf-mean(imsf(:));
imsfdum = imsfdum.*Win;
SF_Pow = fftshift(fftshift(abs(fft2(imsfdum)),1),2);
%figure,imagesc(SF_Pow)

sfradon = radon(imsfdum,[0:10:170]);

[dum idproj] = max(var(odradon));


% odacorr = odradon(:,idproj)/max(odradon(:,idproj));
% sfacorr = sfradon(:,idproj)/max(sfradon(:,idproj));
odacorr = odradon(:,idproj);
sfacorr = sfradon(:,idproj);
% odacorr = xcorr(odacorr,odacorr,'coeff');
% sfacorr = xcorr(sfacorr,sfacorr,'coeff');

dom = (1:length(odacorr))*micperpix;
id = find(odacorr == max(odacorr));
dom = dom-dom(id);

figure
subplot(1,2,1), plot(dom,odacorr/max(odacorr))
hold on,plot(dom,sfacorr/max(sfacorr),'r')
xlim([-1000 1000])

[dum odmin] = min(odacorr);
odmin = dom(odmin);
[dum sfmin] = min(sfacorr);
sfmin = dom(sfmin);
title(['OD min @ ' num2str(odmin) '; SF period @ ' num2str(sfmin) ])
xlabel('microns')

odacorr = odacorr*length(find(imod(:)));
sfacorr = sfacorr*length(find(imod(:)));


odspec_pad = fftshift(abs(fft(odacorr/max(odacorr),65536))); %Pad to interpolate the spectrum (better estimate of peak)
sfspec_pad = fftshift(abs(fft(sfacorr/max(sfacorr),65536)));

dx = dom(2)-dom(1);
fdom_pad = linspace(-1/dx/2,1/dx/2,length(odspec_pad)); %cycles per micron

odspec_pad = odspec_pad(find(fdom_pad>=0));
sfspec_pad = sfspec_pad(find(fdom_pad>=0));
fdom_pad = fdom_pad(find(fdom_pad>=0));

subplot(1,2,2)
plot(fdom_pad,odspec_pad)
hold on,plot(fdom_pad,sfspec_pad,'r')

xlim([0 .003])

[dum id] = max(odspec_pad);
odperiod = 1/fdom_pad(id);
[dum id] = max(sfspec_pad);
sfperiod = 1/fdom_pad(id);

title(['OD period = ' num2str(round(odperiod)) '; SF period = ' num2str(round(sfperiod)) ])
xlabel('cycles per micron')

legend('OD','SF')

%Now recompute spectra w/o padding
odspec = fftshift(abs(fft(odacorr/max(odacorr)))); %Pad to interpolate the spectrum (better estimate of peak)
sfspec = fftshift(abs(fft(sfacorr/max(sfacorr))));
fdom = linspace(-1/dx/2,1/dx/2,length(odspec)); %cycles per micron
odspec = odspec(find(fdom>=0));
sfspec = sfspec(find(fdom>=0));
fdom = fdom(find(fdom>=0));
odspec = odspec(find(fdom>=0));
sfspec = sfspec(find(fdom>=0));
fdom = fdom(find(fdom>=0));

hold on, plot(fdom,odspec,'.b')
hold on,plot(fdom,sfspec,'.r')




%%


% id = find(isnan(imori));
% imori(id) = 0; 
% imoridum = imori-mean(imori(:));
% imoridum = imoridum.*Win;
% ORI_Pow = fftshift(fftshift(abs(fft2(imoridum)),1),2);
% figure,imagesc(ORI_Pow)