function CoMyx = getCellPositions

global maskS cellS ACQinfo


%Resample images to have equal resolution on both axes


masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
celldom = celldom(2:end);

clear CoMyx
for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoMyx(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end
