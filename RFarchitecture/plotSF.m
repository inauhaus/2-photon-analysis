function plotSF(CoM,SF,vthresh)

global ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%10x lens
xmicperpix = xmicperpix*1.5;
ymicperpix = ymicperpix*1.5;

xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

Xmap = SF.pos(2:end);


colorMax = log2(.8);  %Max width on the color bar
colorMin = log2(.01);

idbadX = find(Xmap<colorMin | Xmap>colorMax | SF.varacc(2:end)<vthresh);

Xmap(idbadX) = NaN;

Xmapdum = Xmap;
Xmapdum = (Xmapdum-colorMin)/(colorMax-colorMin);
Xmap_colorIndex = ceil(Xmapdum*64);

cdom = jet;
%%
figure
for p = 1:length(Xmap)
    
    if ~isnan(Xmap(p))
        
        vc = cdom(Xmap_colorIndex(p),:);
        
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end
axis image
colorbar('Ticks',[0 1],'TickLabels',[2^colorMin 2^colorMax])
colormap(cdom)
title('peak SF')
xlabel('microns')

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij

