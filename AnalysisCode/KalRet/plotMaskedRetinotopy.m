function plotMaskedRetinotopy(kmap_hor,kmap_vert) 

%% Plot
%%%%%%%%%%%%%%%
global maskS ACQinfo


sPerx = getparam('x_size');
sPery = getparam('y_size');

%[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
% xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
% ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

xdom = (0:ACQinfo.pixelsPerLine-1);
ydom = (0:ACQinfo.linesPerFrame-1);


masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);

clear CoM
for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

retLim = 140;
cdom = jet;
horplot = (kmap_hor + sPery/2)/sPery; %Should be [0 1] unless there are strange outliers
vertplot = (kmap_vert + sPerx/2)/sPerx;
horplot = horplot*sPery/retLim;
vertplot = vertplot*sPerx/retLim;
horplot = ceil(horplot*64);
vertplot = ceil(vertplot*64);



horplot(find(horplot>64 | horplot<1)) = NaN;
vertplot(find(vertplot>64 | vertplot<1)) = NaN;
figure
for p = 1:length(celldom)
    
    if ~isnan(horplot(p)) & ~isnan(vertplot(p))
       
        if ~isnan(horplot(p)*vertplot(p))
            
            hc = cdom(horplot(p),:);
            vc = cdom(vertplot(p),:);
            
            subplot(2,2,1)
            %plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',hc,'MarkerSize',30)
            %plot(CoM(p,2),CoM(p,1),'o','MarkerFaceColor',hc,'MarkerSize',6,'MarkerEdgeColor','w')
            plot(CoM(p,2),CoM(p,1),'o','MarkerFaceColor',hc,'MarkerSize',6,'MarkerEdgeColor','w')
            hold on
            xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
            
            subplot(2,2,3)
            plot(CoM(p,2),CoM(p,1),'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
            %plot(CoM(p,2)*xmicperpix*ymicperpix,CoM(p,1),'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
            %plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',vc,'MarkerSize',30)
            hold on
            xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
        end
        
    end
end
colormap jet

subplot(2,2,1)
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('horizontal')
colorbar('Ticks',[0 .25 .5 .75 1],'TickLabels',[-retLim/2 -retLim/4 0 retLim/4 retLim/2])
subplot(2,2,3)
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('vertical')
colorbar('Ticks',[0 .25 .5 .75 1],'TickLabels',[-retLim/2 -retLim/4 0 retLim/4 retLim/2])

%% Fit a plane

id = find(isnan(kmap_vert.*kmap_hor));
CoM(id,:) = [];
kmap_hor(id) = [];
kmap_vert(id) = [];

[xmesh ymesh] = meshgrid(xdom,ydom);
H = [CoM(:,2) CoM(:,1) ones(size(CoM,1),1)];

vslope = inv(H'*H)*H'*kmap_vert(:);
vhatIm =  xmesh*vslope(1) + ymesh*vslope(2) + vslope(3);
hslope = inv(H'*H)*H'*kmap_hor(:);
hhatIm =  xmesh*hslope(1) + ymesh*hslope(2) + hslope(3);

subplot(2,2,2),imagesc(xdom,ydom,hhatIm,[-sPery/2 sPery/2]), colorbar, colormap jet, axis image
colorbar('Ticks',[-sPery/2 -sPery/4 0 sPery/4 sPery/2],'TickLabels',[-sPery/2 -sPery/4 0 sPery/4 sPery/2])
hmag = sqrt(sum(hslope(1:2).^2));
title(['Mag = ' num2str(round(1000*hmag)) 'deg/mm'])

subplot(2,2,4),imagesc(xdom,ydom,vhatIm,[-sPery/2 sPery/2]), colorbar, colormap jet, axis image
colorbar('Ticks',[-sPery/2 -sPery/4 0 sPery/4 sPery/2],'TickLabels',[-sPery/2 -sPery/4 0 sPery/4 sPery/2])
vmag = sqrt(sum(vslope(1:2).^2));
title(['Mag = ' num2str(round(1000*vmag)) 'deg/mm'])
xlabel('microns')



%Quantify scatter

vhat =  H*vslope; 
vscat = (vhat-kmap_vert);

hhat =  H*hslope; 
hscat = (hhat-kmap_hor);

R = sqrt(hscat.^2 + vscat.^2);

figure,
hist(R)
xlabel('Distance from plane (deg)')
title(['Median distance from plane = ' num2str(median(R)) 'deg'])






