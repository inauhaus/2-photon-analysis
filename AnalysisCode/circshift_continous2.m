function imshift = circshift_continous2(im,xshift,yshift)

%xshift = 2; yshift = 4;

%Positive x shifts are to the right;
%Positive y shifts are down;

dim = size(im);

sig = .2;

s = round(min(dim)/2);

y = 0:dim(1)-1;
y(end-s+1:end) = -s:1:-1;

x = 0:dim(2)-1;
x(end-s+1:end) = -s:1:-1;



x = x-xshift;

y = y-yshift;



[xdom ydom] = meshgrid(x,y);
r = sqrt(xdom.^2+ydom.^2);
G = exp(-r.^2/(2*sig^2));
G = G/sum(G(:));

imshift = ifft2(fft2(G).*fft2(im));

% figure,imagesc(imshift(100:120,100:120))
% figure,imagesc(im(100:120,100:120))

%figure,imagesc(im)
%figure,imagesc(imshift)