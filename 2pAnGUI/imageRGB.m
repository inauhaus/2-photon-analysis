function imageRGB(im,idx)
%idx is 1,2,or 3, corresponding to R,G,or B
%plot R G B channel

mi = prctile(im(:),.1);
ma = prctile(im(:),99.9);
im(find(im>ma)) = ma;
im(find(im<mi)) = mi;

dim = size(im);
im = (im-min(im(:)));
im = im/max(im(:));

plotter = zeros(dim(1),dim(2),3);
plotter(:,:,idx) = im; 

image((plotter))       %Load and plot frame