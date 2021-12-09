function bwCell = cellMorph(bwCell,mmorph)   

%Morphological opening, followed by downsampling
SE = strel('disk',mmorph,0);
if mmorph ~= 0
    bwCell = imopen(bwCell,SE);
    
   % D = mmorph+1;
    %bwCell = bwCell(1:D:end,1:D:end);
end

