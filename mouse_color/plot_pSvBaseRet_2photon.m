function plot_pSvBaseRet_2photon(pS,kmap_hor,kmap_vert)    

global maskS ACQinfo

%%

Nb = length(pS);
figure,
for bL = 1:Nb
    
    if ~isempty(pS{bL})
        subplot(1,Nb,bL),
        scatter(kmap_vert,pS{bL},'k')
        ylabel('%S')
        xlabel('Vertical retinotopy')
        ylim([0 1])
        xlim([-30 30])
        
        id = find(~isnan(kmap_vert(:).*pS{bL}(:)));
        [r p] = corrcoef(kmap_vert(id),pS{bL}(id));
        r = r(1,2); p = p(1,2);
        
        title(['p = ' num2str(p)])
    end
    
end
    
%% Get Cell locations



[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
celldom = celldom(2:end);

clear CoM
for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end


%% Plot color map
SMmap = jet;
SMmap(:,1) = .4;
SMmap = SMmap(1:40,:);
SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
SMmapdum = SMmap;
SMmap(:,2) = SMmapdum(:,3);
SMmap(:,3) = SMmapdum(:,2);

figure
for bL = 1:Nb
    
    if ~isempty(pS{bL})
        
        pS{bL}(find(pS{bL}>1 | pS{bL}<0)) = NaN;
        pS{bL}(find(isnan(kmap_hor))) = NaN;
        
        pSmap = pS{bL}(2:end);
        pSmap = ceil(pSmap*64);
        cdom = SMmap;
        
        
        subplot(1,Nb,bL)
        for p = 1:length(pSmap)
            
            if ~isnan(pSmap(p))
                
                vc = cdom(pSmap(p),:);
                
                %ax2 = subplot(2,2,4);
                plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',vc,'MarkerSize',30)
                hold on
                xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
            end
        end
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
        title('%Color2 (peak)')
        vR = [0 1];
        %colorbar(ax2,'Ticks',[0 .5 1],'TickLabels',[0 .5 1])
        %colormap(ax2,cdom)
        colorbar('Ticks',[0 .5 1],'TickLabels',[0 .5 1])
        colormap(cdom)
    end
end



