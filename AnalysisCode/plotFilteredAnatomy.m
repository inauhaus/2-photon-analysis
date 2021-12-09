function plotFilteredAnatomy(varargin)

%%
global maskS

%imlog = (maskS.im{1});


%%%%Parameters
rescale = 1;
LPwidth = 0.05;
HPwidth = 100;
black_prc = 1;
white_prc = 99.9
%%%%%%%%%%%%%%%%%

trialdom = [1];
imlog = 0;
for i = 1:length(trialdom)
    dum = GetTrialData([1 0],trialdom(i));
    imlog = imlog + mean(dum{1},3)/length(trialdom);
end

imlog = imlog.^rescale;

h = fspecial('gaussian', size(imlog), LPwidth); h = h/sum(h(:));
imlog = ifft2(abs(fft2(h)).*fft2(imlog));

h = fspecial('gaussian', [size(imlog,1) size(imlog,2)], HPwidth); h = h/sum(h(:));
id = find(h == max(h(:)));
h = -h;
h(id) = h(id)+1;
% h = [h zeros(size(imlog,1),1)];  
% h = [h; zeros(1,size(imlog,1))];
imlog = ifft2(abs(fft2(h)).*fft2(imlog));

figure,
if isempty(varargin)
    imagesc(imlog,[prctile(imlog(:),2) prctile(imlog(:),99.8)]), colormap gray
else
    imagesc(varargin{1},varargin{2},imlog,[prctile(imlog(:),2) prctile(imlog(:),99.8)]), colormap gray
end
axis image

