function imZ = Znormalize(msize)

global maskS

%%%Estimate cell diameter in pixels%%%

%constants
% micpercell = 15; %approximate cell diameter in microns
% micpervolt = 350/2.5; %depends on objective
% 
% micW = micpervolt*ACQinfo.scanAmplitudeX/ACQinfo.zoomFactor;  %FOV micron width 
% micperpixel = micW/ACQinfo.pixelsPerLine; %microns per pixel
%  
% 
% cellD = micperpixel/micpercell;  %Approximate diameter of a cell in pixels
% cellD = cellD*3;



%First truncate outliers
imdum = maskS.im{1};

imdum = imdum - ones(size(imdum,1),1)*mean(imdum);

mi = prctile(imdum(:),.2);
ma = prctile(imdum(:),99.8);
imdum = (imdum-mi)/(ma-mi);
imdum(find(imdum>1)) = 1;
imdum(find(imdum<0)) = 0;

h = fspecial('gaussian',size(imdum),2);
h = h/sum(h(:));
h = abs(fft2(h));


imZ = LocalZ(imdum,round(msize));

imZ = ifft2(fft2(imZ).*h);


figure, imagesc(imZ), colormap gray, colorbar

       
