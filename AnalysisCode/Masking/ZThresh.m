function bw = ZThresh(imZ,Zthresh,mmorph)
    
%It is helpful to interpolate prior to the thesholding and
%morphological opening.  Image is later downsampled by the same amount.
dim = size(imZ);
% if mmorph ~= 0
%     D = 1/(mmorph+1);
%     imZ = interp2(imZ,1:D:dim(2),(1:D:dim(1))');
% end

bw = zeros(size(imZ));
id = find(imZ>Zthresh);
bw(id) = 1;

