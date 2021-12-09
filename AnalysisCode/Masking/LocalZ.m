 function imout = LocalZ(im,W,varargin)
 
%im = im.^2;
 
h = ones(W,W);
[x y] = meshgrid(1:W,1:W);
x = x-(W-1)/2; y = y-(W-1)/2;
rad = sqrt(x.^2 + y.^2);
h(find(rad>(W-1)/2)) = 0; %make it a disc
h = h/sum(h(:));

h = fspecial('gaussian',size(im),W);

mu = conv2(im,h,'same');  %E(x)
mu2 = conv2(im.^2,h,'same');  %E(x^2)
st = sqrt(mu2 - mu.^2); %STD(x)

if isempty(varargin)
    imout = (im-mu)./st;  %Compute local Z value
else
    imout = im-mu;
end

imout = real(imout);