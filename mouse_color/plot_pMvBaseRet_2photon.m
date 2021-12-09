function plot_pMvBaseRet_2photon(pM,kmap_hor,kmap_vert,idExamp)    

global maskS ACQinfo

%%

Nb = length(pM);
figure,
for bL = 1:Nb
    
    if ~isempty(pM{bL})
        subplot(1,Nb,bL),
        scatter(kmap_vert,pM{bL},'k')
        ylabel('%M')
        xlabel('Vertical retinotopy')
        ylim([-.2 1.2])
        xlim([-30 30])
        
        id = find(~isnan(kmap_vert(:).*pM{bL}(:)));
        [r p] = corrcoef(kmap_vert(id),pM{bL}(id));
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
% SMmap = jet;
% SMmap(:,1) = .4;
% SMmap = SMmap(1:40,:);
% SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
% SMmapdum = SMmap;
% SMmap(:,2) = SMmapdum(:,3);
% SMmap(:,3) = SMmapdum(:,2);
% SMmap = flipud(SMmap);

Mmap = zeros(64,3); Mmap(:,2) = (0:63)/63;
Smap = [64 17 64]/64; Smap = (0:63)'*Smap/63;
SMmap = (Smap)+flipud(Mmap);
SMmap = SMmap/max(SMmap(:));
SMmap = flipud(SMmap);

% SMmap = linspace(0,1,64);
% SMmap = SMmap(:)*[1 1 1];
% SMmap(:,2) = flipud(SMmap(:,2));
% SMmap = flipud(SMmap);


figure
for bL = 1:Nb
    
    if ~isempty(pM{bL})
        
        pM{bL}(find(pM{bL}>1.5 | pM{bL}<-.5)) = NaN;
        
        pM{bL}(find(isnan(kmap_hor))) = NaN;
        
        pMmap = pM{bL}(2:end);
        pMmap(find(pMmap>1)) = 1;
        pMmap(find(pMmap<0)) = eps;
        
        pMmap = ceil(pMmap*64);
        cdom = SMmap;
        
        
        subplot(1,Nb,bL)
        for p = 1:length(pMmap)
            
            if ~isnan(pMmap(p))
                
                vc = cdom(pMmap(p),:);
                
                %ax2 = subplot(2,2,4);
                plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','w')
                hold on
                xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image

            end
        end
        
        %Replot the dot (so its on top), and black outline for the examples
        for p = 1:length(pMmap)
            if ~isnan(pMmap(p))
                if ~isempty(find(p == idExamp-1))
                    vc = cdom(pMmap(p),:);
                    
                    %ax2 = subplot(2,2,4);
                    plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'o','MarkerFaceColor',vc,'MarkerSize',6,'MarkerEdgeColor','k')
                    hold on
                    xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image

                end
            end
        end
        
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
        %title('%Color2 (peak)')
        vR = [0 1];
        %colorbar(ax2,'Ticks',[0 .5 1],'TickLabels',[0 .5 1])
        %colormap(ax2,cdom)
        if bL == Nb
            colorbar('Ticks',[0 .5 1],'TickLabels',[0 .5 1])
            colormap(cdom)
        end
        if bL ~= 1
            set(gca,'XTickLabels',[],'YTickLabels',[])
        end
    end
end



