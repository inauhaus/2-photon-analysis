function plotRFpositions(CoM,RF,vthresh)

global ACQinfo

[xmicperpix ymicperpix] = getImResolution;

%10x lens
xmicperpix = xmicperpix*1.5;
ymicperpix = ymicperpix*1.5;

xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

Xmap = RF.xpos(2:end);
Ymap = RF.ypos(2:end);

idbadX = find(Xmap<0 | Xmap>1/getparam('s_freq2') | RF.x_varacc(2:end)<vthresh);
idbadY = find(Ymap<0 | Ymap>1/getparam('s_freq2') | RF.x_varacc(2:end)<vthresh);

Xmap(idbadX) = NaN;
Ymap(idbadY) = NaN;

Xmap_colorIndex = ceil(Xmap*64*getparam('s_freq2'));
Ymap_colorIndex = ceil(Ymap*64*getparam('s_freq2'));

cdom = hsv;
%%
figure
subplot(1,2,1)
for p = 1:length(Xmap)
    
    if ~isnan(Xmap(p))
        
        vc = cdom(Xmap_colorIndex(p),:);
        
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end
axis image
colorbar('Ticks',[0 .5 1],'TickLabels',round([0 1/getparam('s_freq2')/2 1/getparam('s_freq2')]))
colormap(cdom)
title('horizontal retinotopy')
xlabel('microns')

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij

%%

subplot(1,2,2)

for p = 1:length(Ymap)
    
    if ~isnan(Ymap(p))
        
        vc = cdom(Ymap_colorIndex(p),:);
        
        %ax2 = subplot(2,2,4);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
         hold on
         xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        
    end
end

colorbar('Ticks',[0 .5 1],'TickLabels',round([0 1/getparam('s_freq2')/2 1/getparam('s_freq2')]))
colormap(cdom)
title('vertical retinotopy')
xlabel('microns')

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij

axis image
