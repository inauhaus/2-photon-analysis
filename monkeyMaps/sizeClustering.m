function sizeClustering(PW_op,Dbins)

%%
try
    Dist = PW_op.Dist{1}{1};
    
    dori = abs(PW_op.dori{1}{1});
    doripair = PW_op.doripair{1}{1};
    
    dsizepair = PW_op.dsizepair{1}{1};
    dsize = PW_op.dsize{1}{1};

catch
    Dist = PW_op.Dist;
    
    dori = abs(PW_op.dori);   
    doripair = PW_op.doripair;

    
    dsizepair = PW_op.dsizepair;
    dsize = abs(PW_op.dsize);
end
    
clear Distdom dorimu dorisig dsfmu dsfsig SFr SFp dsfBWLinmu dsfBWLinsig SFBWLinr SFBWLinp ORIr OrivSFr OrivSFp Phaser dsizemu dsizesig SIZEr SIZEp 
figure
for i = 1:length(Dbins)-1
    
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dori) & ~isnan(dsize));
   Distdom(i) = mean(Dist(id));
   
   dorimu(i) = median(dori(id));
   dorisig(i) = std(dori(id))/sqrt(length(id));
   
   dsizemu(i) = median(dsize(id));
   dsizesig(i) = std(dsize(id))/sqrt(length(id));  

   [r p] = corrcoef(log2(dsizepair(id,1)),log2(dsizepair(id,2)))
   SIZEr(i) = r(1,2);
   SIZEp(i) = p(1,2);
   
   %[r p] = circCorr(doripair(id,1)*pi/90,doripair(id,2)*pi/90,2*pi,1);
   [r p] = CircCorr_Coherence(doripair(id,1)*pi/90,doripair(id,2)*pi/90,1);
   ORIr(i) = r;
   ORIp(i) = p;
 
Npairs(i) = length(id);
end
Npairs

%%
for i = 1:length(Dbins)-1
    xlab{i} = [num2str(Dbins(i)) '-' num2str(Dbins(i+1))];
end
figure,
subplot(1,2,1)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(ORIr,'-ok'), ylabel('ori1 ori2 correlation')
idSig = find(ORIp<.01);
hold on, plot(idSig,ORIr(idSig),'*r')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)
axis square

subplot(1,2,2)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SIZEr,'-ok'),  ylabel('size size correlation')
idSig = find(SIZEp<.01);
hold on, plot(idSig,SIZEr(idSig),'*r')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)
axis square

%% Plot scatter plot and error bars of pairwise difference
% 
% figure,
% subplot(3,1,1)
% %plot(Dist,dori,'.c')
% scatter(Dist,dori,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
%     'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
% hold on, errorbar(Distdom,dorimu,dorisig,'LineWidth',2), ylabel('dori degrees')
% hold on, plot([0 20],[1 1]*nanmean(dori),'b')
% xlim([0 Dmax])
% title(['r=' num2str(r_oriVDist) '; p=' num2str(p_oriVDist)])
% 
% subplot(3,1,2),
% %scatter(Dist,dsf,'markerEdgeAlpha',.1,'Marker','.','LineWidth',5)
% scatter(Dist,dsf,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
%     'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
% hold on, errorbar(Distdom,dsfmu,dsfsig,'LineWidth',2), ylabel('dsf octaves')
% hold on, plot([0 20],[1 1]*nanmean(dsf),'b')
% xlim([0 Dmax])
% ylim([0 2])
% title(['r=' num2str(r_sfVDist) '; p=' num2str(p_sfVDist)])
% 
% subplot(3,1,3)
% scatter(Dist,dphase,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
%     'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
% hold on, errorbar(Distdom,dphasemu,dphasesig,'LineWidth',2), ylabel('dphase degrees')
% hold on, plot([0 20],[1 1]*nanmean(dphase),'b')
% xlim([0 Dmax])
% title(['r=' num2str(r_phaseVDist) '; p=' num2str(p_phaseVDist)])
% 
% xlabel('mean distance in Bin')
% 
% %%
% 
