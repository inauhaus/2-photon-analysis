 function [imout bwout] = LocalZbw(im,R,thresh)

 global bw
 
dim = size(im);
[rowID colID] = find(double(bw));

[x y] = meshgrid(1:dim(2),1:dim(1));

id = find(bw(:) == 0);
x(id) = NaN;
y(id) = NaN;

imout = NaN*ones(dim(1),dim(2));

for i = 1:length(rowID)
    
    r = sqrt((y-rowID(i)).^2 + (x-colID(i)).^2);
    id = find(r<=R);
    samp = im(id);
    imout(rowID(i),colID(i)) = (imout(rowID(i),colID(i))-mean(samp))/std(samp);
end

bwout = zeros(dim(1),dim(2));
id = find(imout>thresh);
bwout(id) = 1;
