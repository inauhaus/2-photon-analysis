function [bwCell1 bwCell2] = cellMinSize(bwCell,minsize)

%Get rid of cells that are smaller than minsize
for i = 1:2

    if minsize > 0
        
        celllabel = bwlabel(bwCell{i});
        cellid = unique(celllabel);
        for j = 2:length(cellid)
            id = find(cellid(j) == celllabel);
            if length(id) < minsize
                bwCell{i}(id) = 0;
            end
        end
        
    end
    
end

bwCell1 = bwCell{1};
bwCell2 = bwCell{2};