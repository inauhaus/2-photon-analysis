function rst = ColorOri(c1,c2)

%Version 3 uses the function 'analyzeLOGtuning', which generates all of the
%sf parameters.


dfThresh = .05;

global cellS idExamp ACQinfo Analyzer

%Order in the looper

for p = 1:length(Analyzer.L.param)
    switch Analyzer.L.param{p}{1}
        case 'theta'
            colordim = p;
        case 'ori'
            oridim = p;            
    end
end

[dum dimsort] = sort([colordim oridim]);


%ORI domains
oridom = getdomain('ori');

% oridom = oridom(1:end/2);
% oridomI = linspace(oridom(1),180,length(oridom)*2); %oridomI = oridomI(1:end-1);
% oridomII = linspace(oridom(1),180,length(oridom)*5); %oridomII = oridomII(1:end-1);%strictly for plotting a smooth curve

%These are used for plotting outside of this function
rst.oridom = oridom;
% rst.oridomI = oridomI;
% rst.oridomII = oridomII;

cid = [c1 c2];

colordom = getdomain('theta');
%figure
for i = 1:length(cellS.mukern)

    %kdum = abs(cellS.F1kern{i});
    kdum = cellS.mukern{i};
    kdum = permute(kdum,dimsort);  %force the kernel to have dimensions in this order
    
    kdumsig = cellS.sigkern{i};
    kdumsig = permute(kdumsig,dimsort);  %force the kernel to have dimensions in this orde
    
    kdumT = squeeze(cellS.mukernTime{i});
    kdumT = permute(kdumT,[dimsort 3]);
    
    %sforimat = squeeze(mean(kdum,1)); %mean over color
    
    oricolormat = kdum;
    rst.oricolormat(:,:,i) = oricolormat;
    
    oricolormat_sig = kdumsig;
    rst.oricolormat_sig(:,:,i) = oricolormat_sig; %Used for example plots
    
    rst.kernT(:,:,:,i) = kdumT;  %Used for example plots
    
    %Compute tuning curves using a weighted average of the orthogonal
    %dimension
%     cdum = sum(oricolormat);
%     cdum  = cdum/sum(cdum);
%     oritc = sum(sforimat.*(ones(length(oridom),1)*cdum),2);
    oritc = mean(oricolormat,1);
    
%     oridum = sum(oricolormat,2);
%     oridum = oridum/sum(oridum);
%     colortc = sum(oricolormat.*(oridum*ones(1,length(colordom))),1);
    colortc = mean(oricolormat,2);
    
    
 
    
    %pick a condition for the example time course
    %         sforimatdum = squeeze(kdum(2,:,:)); %use one color for ori and sf
    %         [idy idx] = find(sforimatdum == max(sforimatdum(:)));
    %
    %         rst.tcourse{c}(i,:) = squeeze(kdumT(c,idy,idx,:));
    
    %         idy2 = idy+length(oridom)/2;
    %         if idy2>length(oridom)
    %             idy2 = idy2-length(oridom);
    %         end
    %         rst.tcourse{c}(i,:) = (rst.tcourse{c}(i,:) + squeeze(kdumT(c,idy2,idx,:))')/2;

    
%     subplot(ceil(sqrt(length(cellS.mukern))),ceil(sqrt(length(cellS.mukern))),i)
%     %         hold on
%     plot(oridom,oricolormat(1,:),'.-g')
%     hold on
%     plot(oridom,oricolormat(2,:),'.-b')
%     ylim([0 .2])
    
    
    %hold on
    %plot(oridom,oritc)
    
    
    %         plot(rst.tcourse{c}(i,:))
    %         hold on
    
    %title(num2str(sfpref(i)))
     
    %colortc = phi(colortc);
    
    idM = find(colordom == 0);
    idS = find(colordom == 90);
    rst.pM(i) = colortc(idM)/(colortc(idM)+colortc(idS));
    
    rst.M_Resp(i) = colortc(idM);
    rst.S_Resp(i) = colortc(idS);
    
%      if rst.pS(i) < 0
%          rst.pS(i) = NaN;
%      end
    
%rst.pS(i) = phi(rst.pS(i));


% kdumT = squeeze(mean(kdumT,2));
% tc1 = kdumT(1,:);
% tc2 = kdumT(3,:);

   if max(cellS.mukern{i}(:))<dfThresh
        
       rst.kernT(:,:,:,i) = rst.kernT(:,:,:,i)*NaN;
       rst.oricolormat_sig(:,:,i) = rst.oricolormat_sig(:,:,i)*NaN;
       rst.oricolormat(:,:,i) = rst.oricolormat(:,:,i)*NaN;
       rst.pM(i) = rst.pM(i)*NaN;
       
       rst.M_Resp(i) = rst.M_Resp(i)*NaN;
       rst.S_Resp(i) = rst.S_Resp(i)*NaN;
       
    end

    
    
end

%%

%oridomplot = [oridom 180];

% 
% fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
% tdom =(0:(size(rst.tcourse{1},2)-1))*fr/1000; 
% figure
% for i = 1:length(idExamp)
%     
%     subplot(length(idExamp),5,(i-1)*5+1)
%     imagesc(sfdom,oridom,rst.sforimat{1}(:,:,idExamp(i))), axis square
%     set(gca,'YTick',[0 180])
%     if i == length(idExamp)
%         xlabel('SF (cyc/deg)')
%         ylabel('Direction (deg)')
%     end
%     if i == 1
%        title('color 1') 
%     end
%     
%     subplot(length(idExamp),5,(i-1)*5+2)
%     imagesc(sfdom,oridom,rst.sforimat{2}(:,:,idExamp(i))), axis square
%     set(gca,'YTick',[0 180])
%     colormap gray    
%     if i == 1
%        title('color 2') 
%     end
%      
%     subplot(length(idExamp),5,(i-1)*5+3)
%     
%     polarplot([oridom 360]*pi/180,[rst.oritc{1}(idExamp(i),:) rst.oritc{1}(idExamp(i),1)],'.-k')
%     hold on
%     polarplot([oridom 360]*pi/180,[rst.oritc{2}(idExamp(i),:) rst.oritc{2}(idExamp(i),1)],'.-r')
%     hold on,
%     polarplot([0 0],[0 0],'.k')
%     axis off
%     %ax = gca;
%     %ax.ThetaTick = [0 90 180 270];
%     %set(gca,'ThetaTick',[0 90 180 270])
%     %hold on
%     %plot([oridomII 180],[rst.oritcfit{2}(idExamp(i),:) rst.oritcfit{2}(idExamp(i),1)],'r')  
% %     xlim([0 360])
% %     axis tight
% %     yl = ylim;
% %     yl(1) = min([0 yl(1)]);
% %     ylim(yl)
% %     set(gca,'XTick',[0 90 180])
%     
%     
% %     
%     subplot(length(idExamp),5,(i-1)*5+5)
%     plot(tdom,rst.tcourse{1}(idExamp(i),:),'k')
%     hold on
%     plot(tdom,rst.tcourse{2}(idExamp(i),:),'r')
%     hold on
%     plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
%     xlabel('sec'), ylabel('dF/F')
%     axis tight
%     %title(num2str(sfpref(i)))
% 
%     
% end

%%
% 
% global maskS 
% 
% [xmicperpix ymicperpix] = getImResolution;
% 
% %Resample images to have equal resolution on both axes
% xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
% ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
% 
% masklabel = bwlabel(maskS.bwCell{1},4);
% celldom = unique(masklabel);
% celldom = celldom(2:end);
% 
% clear CoM
% for p = 1:length(celldom)
%     [idcelly idcellx] = find(masklabel == celldom(p));
%     CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
% end
% 
% %% Plot color map
% SMmap = jet;
% SMmap(:,1) = .4;
% SMmap = SMmap(1:40,:);
% SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
% SMmapdum = SMmap;
% SMmap(:,2) = SMmapdum(:,3);
% SMmap(:,3) = SMmapdum(:,2);
% 
% 
% rst.pS(find(rst.pS>1 | rst.pS<0)) = NaN;
% 
% pSmap = rst.pS(2:end);
% pSmap = ceil(pSmap*64);
% cdom = SMmap;
% 
% figure
% for p = 1:length(pSmap)
%     
%     if ~isnan(pSmap(p))
%         
%         vc = cdom(pSmap(p),:);
%         
%         %ax2 = subplot(2,2,4);
%         plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',vc,'MarkerSize',30)
%         hold on
%         xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
%     end
%     
% end
% 
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
% title('%Color2 (peak)')
% vR = [0 1];
% %colorbar(ax2,'Ticks',[0 .5 1],'TickLabels',[0 .5 1])
% %colormap(ax2,cdom)
% colorbar('Ticks',[0 .5 1],'TickLabels',[0 .5 1])
% colormap(cdom)
% 
