function manualAddCell_GCaMP

global fh maskS

fh = gcf;

[xmicperpix ymicperpix] = getImResolution(0);
res = geomean([xmicperpix ymicperpix]);

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);

%Bandpass filter
im = maskS.im{1};

LPsig = 1*(1.6/res);
HPsig = 2*(1.6/res);

LP = fspecial('gaussian',size(im),1); LP = LP/sum(LP(:));
HP = fspecial('gaussian',size(im),2); HP = HP/sum(HP(:));
BP = LP-HP;
imH = ifft2(fft2(BP).*fft2(im));
imH = fftshift(fftshift(imH,1),2);
maskS.BP = imH;


