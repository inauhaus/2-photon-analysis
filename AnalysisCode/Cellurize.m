function [im locs] = Cellurize(im,bwCell)

global bw locx locy

dim = size(im);

cellid = bwlabel(bwCell.*bw);
uid = unique(cellid);
uid = uid(2:end); %First element is the background (0)

for i = 1:length(uid)  
    idx = find(cellid(:) == uid(i));
    im(idx) = mean(im(idx));

    locxdum = ceil(idx/dim(1));
    locydum = idx-dim(1)*(locxdum-1);
    locx(i) = round(mean(locxdum));
    locy(i) = round(mean(locydum));
    
end

locs = (locx-1)*dim(1) + locy;

%     id = find(~bwCell);
%     im(id) = 0;
