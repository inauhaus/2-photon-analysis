function plotDotMap(vals,cdom,masklabel,celldom,valRange)

global maskS ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

clear CoM
for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

vals = (vals-valRange(1))/(valRange(2)-valRange(1));

vals = round(vals*63 +1);

vals(find(vals>64)) = 64;
vals(find(vals<1)) = 1;

if isstr(cdom)  %e.g. if cdom = 'jet'
    cdom = eval(cdom);
end

for p = 1:length(celldom)
    
    if ~isnan(vals(p))
        
        RGB = cdom(vals(p),:);
        
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',RGB,'MarkerSize',30)
        hold on
              
    end
end

axis image
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), 
axis ij
