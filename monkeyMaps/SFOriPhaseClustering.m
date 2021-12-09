function SFOriPhaseClustering(PW_ro,Dbins)

%%
try
    Dist = PW_ro.Dist{1};
    dori = abs(PW_ro.dori{1});
    dsf = abs(PW_ro.dsf{1});
    dsfBWLin = abs(PW_ro.dsfBWLin{1});
    dphase = abs(PW_ro.dphaseA{1});
    dF1F0 = abs(PW_ro.dF1F0{1});
    dsfpair = PW_ro.dsfpair{1};
    dsfBWLinpair = PW_ro.dsfBWLinpair{1};  %sf BW (2sigma)
    dsfBWpair = ((dsfBWLinpair/2 + dsfpair)./dsfpair);
    doripair = PW_ro.doripair{1};
    dorisigpair = PW_ro.dorisigpair{1};
    dphasepair = PW_ro.dphaseApair{1};
    dF1F0pair = PW_ro.dF1F0pair{1};
    dLPnesspair = PW_ro.dLPnesspair{1};

catch
    Dist = PW_ro.Dist;
    dori = abs(PW_ro.dori);
    dsf = abs(PW_ro.dsf);
    dsfBWLin = abs(PW_ro.dsfBWLin);
    dphase = abs(PW_ro.dphaseA);
    dF1F0 = abs(PW_ro.dF1F0);
    dsfpair = PW_ro.dsfpair;
    dsfBWLinpair = PW_ro.dsfBWLinpair;
    dsfBWpair = ((PW_ro.dsfBWLinpair/2 + PW_ro.dsfpair)./PW_ro.dsfpair);
    doripair = PW_ro.doripair;
    dorisigpair = PW_ro.dorisigpair;
    dphasepair = PW_ro.dphaseApair;
    dF1F0pair = PW_ro.dF1F0pair;
    dLPnesspair = PW_ro.dLPnesspair;
end
    

doriBW = abs(dorisigpair(:,1)-dorisigpair(:,2));
dsfBW = abs(log2(dsfBWpair(:,1)./dsfBWpair(:,2)));

%Dbins = 0:20:100;
clear Distdom Npairs
clear dsfmu dsfsig SFr SFp dsfBWLinmu dsfBWLinsig SFBWLinr SFBWLinp 
clear dsfBWLinmu dsfBWLinsig SFBWLinr SFBWLinp
clear dsfBWmu dsfBWsig SFBWr SFBWp
clear dorimu dorisig doriBWmu doriBWsig ORIr OrivSFr OrivSFp 
clear Phaser Phaserp dphasemu dphasesig

figure
for i = 1:length(Dbins)-1
    
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dori) & ~isnan(dsf));
   Distdom(i) = mean(Dist(id));
   
   dorimu(i) = median(dori(id));
   dorisig(i) = std(dori(id))/sqrt(length(id));
   
   doriBWmu(i) = median(doriBW(id));
   doriBWsig(i) = std(doriBW(id))/sqrt(length(id));
   
   dsfmu(i) = median(dsf(id));
   dsfsig(i) = std(dsf(id))/sqrt(length(id));
   
   dsfBWLinmu(i) = median(dsfBWLin(id));
   dsfBWLinsig(i) = std(dsfBWLin(id))/sqrt(length(id));
   
   dsfBWmu(i) = median(dsfBW(id));
   dsfBWsig(i) = std(dsfBW(id))/sqrt(length(id));
   
   dphasemu(i) = median(dphase(id));
   dphasesig(i) = std(dphase(id))/sqrt(length(id));
   
   [r p] = corrcoef(log2(dsfBWLinpair(id,1)),log2(dsfBWLinpair(id,2)));
   SFBWLinr(i) = r(1,2);
   SFBWLinp(i) = p(1,2);
   
   [r p] = corrcoef(log2(dsfBWpair(id,1)),log2(dsfBWpair(id,2)));
   SFBWr(i) = r(1,2);
   SFBWp(i) = p(1,2);
   
   Npairs(i) = length(id)
   [r p] = corrcoef(log2(dsfpair(id,1)),log2(dsfpair(id,2)));
   SFr(i) = r(1,2);
   SFp(i) = p(1,2);
   
   [r p] = corrcoef(log2(dorisigpair(id,1)),log2(dorisigpair(id,2)));
   oriBWr(i) = r(1,2);
   oriBWp(i) = p(1,2);
   
   %[r p] = circCorr(doripair(id,1)*pi/90,doripair(id,2)*pi/90,2*pi,1);
   [Proj p] = CircCorr_Coherence(doripair(id,1)*pi/90,doripair(id,2)*pi/90,1);
   ORIr(i) = Proj;
   ORIp(i) = p;
   
   %[r p] = circCorr(dphasepair(id,1)*pi/90,dphasepair(id,2)*pi/180,2*pi,1);
   [r p] = CircCorr_Coherence(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180,1);
   Phaser(i) = r;
   Phaserp(i) = p;
   
   idF1F0nan = find(~isnan(dF1F0pair(id,1).*dF1F0pair(id,2)) & dLPnesspair(id,1)<0.61 & dLPnesspair(id,2)<0.61);
   [r p] = corrcoef((dF1F0pair(id(idF1F0nan),1)),(dF1F0pair(id(idF1F0nan),2)));
   F1F0r(i) = r(1,2);
   F1F0p(i) = p(1,2);
   
   [r p] = corrcoef(dori(id),dsf(id));
   OrivSFr(i) = r(1,2);
   OrivSFp(i) = p(1,2);
   
   subplot(1, length(Dbins),i), 
   scatter(dori(id),dsf(id),'.')
   [r p] = corrcoef(dori(id),dsf(id));
   title(['r=' num2str(r(1,2)) 'p=' num2str(p(1,2))])
   
   
   
   doriBinEdges = [0:22.5:90];
   clear mudori mudsf stdori stdsf
   for i = 1:length(doriBinEdges)-1
       
       id2 = find(dori(id)>doriBinEdges(i) & dori(id)<doriBinEdges(i+1));
       
       mudori(i) = median(dori(id(id2)));
       mudsf(i) = mean(dsf(id(id2)));
       
       stddori(i) = std(dori(id(id2)));
       stddsf(i) = std(dsf(id(id2)));
              
   end
   
   hold on, errorbar(mudori,mudsf,stddsf,'r');
   xlabel('dori'), ylabel('dsf')
   
   
end

%%

%phase will have fewer pairs because they will often have different
%preferred ori/sf, so use a different conditional for pairings...
clear dphasemu dphasesig Phaser Distdom_phase
%BPpair = find(dLPnesspair(:,1)>.5 & dLPnesspair(:,2)>.5);  %Only use pairs that are more bandpass
for i = 1:length(Dbins)-1
    
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dphase));
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dphase) & dLPnesspair(:,1)<.5 & dLPnesspair(:,2)<.5);
   Distdom_phase(i) = mean(Dist(id));
   
   dphasemu(i) = median(dphase(id));
   dphasesig(i) = std(dphase(id))/sqrt(length(id));
   
   [r p] = circCorr(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180,2*pi,1);
   [r p] = CircCorr_Coherence(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180,1);
   Phaser(i) = r;
   Phaserp(i) = p;
   
end

%% Phase clustering for different preferred sf

for i = 1:length(Dbins)-1
    xlab{i} = [num2str(Dbins(i)) '-' num2str(Dbins(i+1))];
end

sfBin = [0 1 2 4];
musf = sqrt(dsfpair(:,1).*dsfpair(:,2));

%phase will have fewer pairs because they will often have different
%preferred ori/sf, so use a different conditional for pairings...
clear dphasemu2 dphasesig2 Phaser_band Phaserp_band Distdom_phase2
for i = 1:length(Dbins)-1
    for j = 1:length(sfBin)-1
        

        id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dphase) & musf>sfBin(j) & musf<sfBin(j+1));
        Distdom_phase2(i,j) = mean(Dist(id));
        
        dphasemu2(i,j) = median(dphase(id));
        dphasesig2(i,j) = std(dphase(id))/sqrt(length(id));
        
        %[r p] = circCorr(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180,2*pi,1);
        [r p] = CircCorr_Coherence(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180,1);
        Phaser_band(i,j) = r;
        Phaserp_band(i,j) = p;
        
        
    end
end

figure, plot(Phaser_band)
set(gca,'XTickLabel',xlab,'FontSize',7)
legend('lowSF','midSF','hiSF')

%% Get correlation coeffs at close range


id = find(Dist>0 & Dist<Dbins(end) & ~isnan(dori) & ~isnan(dsf));
[r p] = corrcoef(Dist(id),dori(id));
r_oriVDist = r(1,2);
p_oriVDist = p(1,2);

[r p] = corrcoef(Dist(id),dsf(id));
r_sfVDist = r(1,2);
p_sfVDist = p(1,2);

id = find(Dist>0 & Dist<Dbins(end) & ~isnan(dphase));
[r p] = corrcoef(Dist(id),dphase(id));
r_phaseVDist = r(1,2);
p_phaseVDist = p(1,2);


%%
for i = 1:length(Dbins)-1
    xlab{i} = [num2str(Dbins(i)) '-' num2str(Dbins(i+1))];
end
figure,
subplot(2,4,1)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(ORIr,'-ok'), ylabel('ori1 ori2 correlation')
idSig = find(ORIp<.01);
hold on, plot(idSig,ORIr(idSig),'*r')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)
axis square

subplot(2,4,2)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFr,'-ok'),  ylabel('sf1 sf2 correlation')
idSig = find(SFp<.01);
hold on, plot(idSig,SFr(idSig),'*r')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)
axis square

subplot(2,4,3)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFBWLinr,'-ok'),  ylabel('sfBW1 sfBW2 correlation')
idSig = find(SFBWLinp<.01);
hold on, plot(idSig,SFBWLinr(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
axis square

subplot(2,4,4)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFBWr,'-ok'),  ylabel('sflogBW1 sflogBW2 correlation')
idSig = find(SFBWp<.01);
hold on, plot(idSig,SFBWr(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
axis square

subplot(2,4,5)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(OrivSFr,'-ok'), ylabel('|dori| vs. |dsf| correlation')
idSig = find(OrivSFp<.01);
hold on, plot(idSig,OrivSFr(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
axis square
ylim([-.5 1])

subplot(2,4,6)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(Phaser,'-ok'), ylabel('phase1 phase2 correlation')
idSig = find(Phaserp<.01);
hold on, plot(idSig,Phaser(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
xlabel('distance bins (um)')
axis square
ylim([-.5 1])

subplot(2,4,7)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(oriBWr,'-ok'),  ylabel('oriBW1 oriBW2 correlation')
idSig = find(oriBWp<.01);
hold on, plot(idSig,oriBWr(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
axis square

subplot(2,4,8)
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(F1F0r,'-ok'),  ylabel('F1F01 F1F02 correlation')
idSig = find(F1F0p<.01);
hold on, plot(idSig,F1F0r(idSig),'*r')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
axis square

%%
figure
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFr,'-ok','Color',[0 0 1]),  ylabel('sf1 sf2 correlation')
idSig = find(SFp<.01);
hold on, plot(idSig,SFr(idSig),'*b')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)

plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFBWLinr,'-ok','Color',[0 0 0]),  ylabel('sfBW1 sfBW2 correlation')
idSig = find(SFBWLinp<.01);
hold on, plot(idSig,SFBWLinr(idSig),'*k')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])

plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFBWr,'--ok','Color',[0 0 0]),  ylabel('sfBW1 sfBW2 correlation')
idSig = find(SFBWp<.01);
hold on, plot(idSig,SFBWr(idSig),'*k')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])

figure
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFr,'-ok','Color',[0 0 1]),  ylabel('sf1 sf2 correlation')
idSig = find(SFp<.01);
hold on, plot(idSig,SFr(idSig),'*b')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)

plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(oriBWr,'-ok','Color',[0 0 0]),  ylabel('oriBW1 oriBW2 correlation')
idSig = find(oriBWp<.01);
hold on, plot(idSig,oriBWr(idSig),'*k')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
title('oriBW vs SF')

figure
plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(SFr,'-ok','Color',[0 0 1]),  ylabel('sf1 sf2 correlation')
idSig = find(SFp<.01);
hold on, plot(idSig,SFr(idSig),'*b')
ylim([-.5 1])
set(gca,'XTickLabel',xlab,'FontSize',7)

plot([1 length(Dbins)-1],[0 0],'--','Color',[.6 .6 .6])
hold on
plot(F1F0r,'-ok','Color',[0 0 0]),  ylabel('F1F01 F1F02 correlation')
idSig = find(F1F0p<.01);
hold on, plot(idSig,F1F0r(idSig),'*k')
set(gca,'XTickLabel',xlab,'FontSize',7)
ylim([-.5 1])
title('F1F0 vs SF')
%% Plot scatter plot and error bars of pairwise difference



figure,
subplot(3,1,1)
%plot(Dist,dori,'.c')
scatter(Dist,dori,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
hold on, errorbar(Distdom,dorimu,dorisig,'LineWidth',2), ylabel('dori degrees')
hold on, plot([0 20],[1 1]*nanmean(dori),'b')
xlim([0 Dbins(end)])
title(['r=' num2str(r_oriVDist) '; p=' num2str(p_oriVDist)])

subplot(3,1,2),
%scatter(Dist,dsf,'markerEdgeAlpha',.1,'Marker','.','LineWidth',5)
scatter(Dist,dsf,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
hold on, errorbar(Distdom,dsfmu,dsfsig,'LineWidth',2), ylabel('dsf octaves')
hold on, plot([0 20],[1 1]*nanmean(dsf),'b')
xlim([0 Dbins(end)])
ylim([0 2])
title(['r=' num2str(r_sfVDist) '; p=' num2str(p_sfVDist)])

subplot(3,1,3)
scatter(Dist,dphase,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
hold on, errorbar(Distdom,dphasemu,dphasesig,'LineWidth',2), ylabel('dphase degrees')
hold on, plot([0 20],[1 1]*nanmean(dphase),'b')
xlim([0 Dbins(end)])
title(['r=' num2str(r_phaseVDist) '; p=' num2str(p_phaseVDist)])

xlabel('mean distance in Bin')

%%

dsfLin = abs(dsfpair(:,1)-dsfpair(:,2));
varsf = var(dsfpair,[],2);
clear mudsf Distdom
for i = 1:length(Dbins)-1
    
    
    %id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dsfLin));
    id = find(Dist<Dbins(i+1) & ~isnan(dsfLin));
    Distdom(i) = mean(Dist(id));
    mudsf(i) = mean(varsf(id)).^.5;
    
    
end

figure,plot(Distdom,mudsf,'-o')
xlabel('Distance'), ylabel('Estimate of std in pooling window')

