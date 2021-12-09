function plotArchitecture(CoM,RF,vthresh)

global ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%10x lens
xmicperpix = xmicperpix*1.5;
ymicperpix = ymicperpix*1.5;

xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

Xmap = RF.xsig(2:end);
Ymap = RF.ysig(2:end);
Amap = sqrt(RF.ysig(2:end).^2 + RF.xsig(2:end).^2);

colorMax = 25;  %Max width on the color bar

idbadX = find(Xmap<2 | Xmap>colorMax | RF.x_varacc(2:end)<vthresh);
idbadY = find(Ymap<2 | Ymap>colorMax | RF.y_varacc(2:end)<vthresh);

Xmap(idbadX) = NaN;
Ymap(idbadY) = NaN;
Amap(unique([idbadX(:); idbadY(:)])) = NaN;

Xmap_colorIndex = ceil(Xmap*64/colorMax);
Ymap_colorIndex = ceil(Ymap*64/colorMax);
Amap_colorIndex = ceil(Amap*64/colorMax);

Amap

cdom = jet;
%%
figure
subplot(1,3,1)
for p = 1:length(Xmap)
    
    if ~isnan(Xmap(p))
        
        vc = cdom(Xmap_colorIndex(p),:);
        
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end
axis image
colorbar('Ticks',[0 .5 1],'TickLabels',[0 colorMax/2 colorMax])
colormap(cdom)
title('horizontal width (sigma)')
xlabel('microns')

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij

%%

subplot(1,3,2)

for p = 1:length(Ymap)
    
    if ~isnan(Ymap(p))
        
        vc = cdom(Ymap_colorIndex(p),:);
        
        %ax2 = subplot(2,2,4);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
         hold on
         xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end

title('vertical width (sigma)')
xlabel('microns')

colorbar('Ticks',[0 .5 1],'TickLabels',[0 colorMax/2 colorMax])
colormap(cdom)

axis image

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij

%%

subplot(1,3,3)

for p = 1:length(Amap)
    
    if ~isnan(Amap(p))
        
        vc = cdom(Amap_colorIndex(p),:);
        
        %ax2 = subplot(2,2,4);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
         hold on
         xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end

title('sqrt(Area) (sigma)')
xlabel('microns')

colorbar('Ticks',[0 .5 1],'TickLabels',[0 colorMax/2 colorMax])
colormap(cdom)

axis image

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij