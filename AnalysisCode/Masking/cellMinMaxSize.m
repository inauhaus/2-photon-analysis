function bwCell = cellMinMaxSize(bwCell,mima)


if mima(1) > 0
    
    celllabel = bwlabel(bwCell,4);  % 4 so that corners are not considered "connected"
    cellid = unique(celllabel);
    for j = 1:length(cellid)
        id = find(cellid(j) == celllabel);
        if length(id) < mima(1) || length(id) > mima(2)
            bwCell(id) = 0;
        end
    end
    
end
