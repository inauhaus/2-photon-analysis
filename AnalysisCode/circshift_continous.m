function imshift = circshift_continous(im,xshift,yshift)

%xshift = 2; yshift = 4;

%Positive x shifts are to the right;
%Positive y shifts are down;

PadFactor = 2;

dim = size(im);
dimPad = dim*PadFactor;

% wdomx = linspace(0,2*pi,size(im,2)+1);
% wdomy = linspace(0,2*pi,size(im,1)+1);
% wdomx = wdomx(1:end-1);
% wdomy = wdomy(1:end-1);

wdomx = linspace(0,2*pi,dim(2));
wdomy = linspace(0,2*pi,dim(1));
wdomx = wdomx(1:end);
wdomy = wdomy(1:end);

[wx wy] = meshgrid(wdomx,wdomy);

imF = fft2(im);
%imF = imF.*exp(-1i*wx*xshift).*exp(-1i*wy*yshift);
imF = imF.*exp(-1i*wx*(xshift));

imFPad = repmat(imF,2);

imshift = real(ifft2(imFPad));

% figure,imagesc(imshift)
% figure,imagesc(im)

figure,imagesc(imshift(100:120,100:120))
figure,imagesc(im(100:120,100:120))