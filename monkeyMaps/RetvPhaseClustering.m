function RetvPhaseClustering(PW_ro,PW_op)

%%
try
    Dist = PW_ro.Dist{1};
    dori = abs(PW_ro.dori{1});
    dsf = abs(PW_ro.dsf{1});
    dphase = abs(PW_ro.dphaseA{1});
    dsfpair = PW_ro.dsfpair{1};
    doripair = PW_ro.doripair{1};
    dphasepair = PW_ro.dphaseApair{1};
    
    muori = exp(1i*PW_ro.doripair{1}(:,1)*pi/180*2) + exp(1i*PW_ro.doripair{1}(:,2)*pi/180*2);
    muori = angle(muori)*pi/180/2;
    muori(find(muori<0)) = muori(find(muori<0))+180;
    dpos = PW_op.dpos{1}{1};
    ax = PW_op.ax{1}{1};
    dx = dpos.*cos(ax*pi/180);
    dy = dpos.*sin(ax*pi/180);
    dposOridim = abs(sum([dx dy].*[cos(muori*pi/180) sin(muori*pi/180)],2));
    
catch
    Dist = PW_ro.Dist;
    dori = abs(PW_ro.dori);
    dsf = abs(PW_ro.dsf);
    dphase = (PW_ro.dphaseA);
    dsfpair = PW_ro.dsfpair;
    doripair = PW_ro.doripair;
    dphasepair = PW_ro.dphaseApair;
    
    muori = exp(1i*PW_ro.doripair(:,1)*pi/180*2) + exp(1i*PW_ro.doripair(:,2)*pi/180*2);
    muori = angle(muori)*pi/180/2;
    muori(find(muori<0)) = muori(find(muori<0))+180;
    dpos = PW_op.dpos;
    ax = PW_op.ax;
    dx = dpos.*cos(ax*pi/180);
    dy = dpos.*sin(ax*pi/180);
    dposOridim = abs(sum([dx dy].*[cos(muori*pi/180) sin(muori*pi/180)],2));
end

musf = sqrt(dsfpair(:,1).*dsfpair(:,2));
id = find(musf>2);
dphase(id) = NaN;

degphaseshift = dphase(:)./musf(:)/360;



%%
maxD = 200;


figure,
id = find(~isnan(Dist.*dphase) & Dist<maxD)
scatter(Dist(id),dphase(id),'.')
[r p] = corrcoef(Dist(id),dphase(id))
xlabel('dpos (mm of cortex)')
ylabel('d phase')

figure,
id = find(~isnan(dpos.*dphase) & Dist<maxD)
scatter(dpos(id),dphase(id)./musf(id)/360,'.')
[r p] = corrcoef(dpos(id),dphase(id)./musf(id))
xlabel('d retinotopy')
ylabel('d phase/360/dori')


figure,
id = find(~isnan(Dist.*dpos) & Dist<maxD)
scatter(Dist(id),dpos(id),'.')
[r p] = corrcoef(Dist(id),dpos(id))
xlabel('dpos (mm of cortex)')
ylabel('d retinotopy')

%%

maxD = 200;

figure,

id = find(~isnan(dposOridim.*degphaseshift) & Dist<maxD)
subplot(2,2,2)
scatter(dposOridim(id),degphaseshift(id),'.k');
axis square
xlim([0 .6])
ylim([0 .6])
ylabel('phase shift / sfreq (deg)')
xlabel('RF shift (deg)')
[r p] = corrcoef(dposOridim(id),degphaseshift(id))
mu = round(median(degphaseshift(id))*100)/100;
subplot(2,2,1)
hist(degphaseshift(id),[0:.1:.8]), xlabel('phase shift / sfreq (deg)')
title(['mean = ' num2str(mu) ' deg'])
ylabel('n pairs')
subplot(2,2,4)
hist(dposOridim(id),[0:.1:.8])
mu = round(median(dposOridim(id))*100)/100;
title(['mean = ' num2str(mu) ' deg'])
xlabel('RF shift (deg)')
ylabel('n pairs')

%%
Dmax = 400;
%dsf = abs(log2(PW_ro.dsfpair{1}(:,1)./PW_ro.dsfpair{1}(:,2)));

Dbins = 0:50:300;
%Dbins = 0:20:100;
clear Distdom dorimu dorisig dsfmu dsfsig SFr SFp ORIr OrivSFr OrivSFp Phaser
figure
for i = 1:length(Dbins)-1
    
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dori) & ~isnan(dsf));
   Distdom(i) = mean(Dist(id));
   
   dorimu(i) = mean(dori(id));
   dorisig(i) = std(dori(id))/sqrt(length(id));
   
   dsfmu(i) = mean(dsf(id));
   dsfsig(i) = std(dsf(id))/sqrt(length(id));
   
   dphasemu(i) = mean(dphase(id));
   dphasesig(i) = std(dphase(id))/sqrt(length(id));
   
   [r p] = corrcoef(log2(dsfpair(id,1)),log2(dsfpair(id,2)));
   SFr(i) = r(1,2);
   SFp(i) = p(1,2);
   
   [r] = circCorr(doripair(id,1)*pi/90,doripair(id,2)*pi/90);
   ORIr(i) = r;
   %ORIp(i) = p(1,2);
   
   [r] = circCorr(dphasepair(id,1)*pi/90,dphasepair(id,2)*pi/90);
   Phaser(i) = r;
   %ORIp(i) = p(1,2);
   
   [r p] = corrcoef(dori(id),dsf(id));
   OrivSFr(i) = r(1,2);
   OrivSFp(i) = p(1,2);
   
   subplot(1, length(Dbins),i), 
   scatter(dori(id),dsf(id),'.')
   [r p] = corrcoef(dori(id),dsf(id))
   title(['r=' num2str(r(1,2)) 'p=' num2str(p(1,2))])
end

%phase will have fewer pairs because they will often have different
%preferred ori/sf, so use a different conditional for pairings...
clear dphasemu dphasesig Phaser Distdom_phase
for i = 1:length(Dbins)-1
    
   id = find(Dist>Dbins(i) & Dist<Dbins(i+1) & ~isnan(dphase));
   Distdom_phase(i) = mean(Dist(id));
   
   dphasemu(i) = median(dphase(id));
   dphasesig(i) = std(dphase(id))/sqrt(length(id));
   
   [r] = circCorr(dphasepair(id,1)*pi/180,dphasepair(id,2)*pi/180);
   Phaser(i) = r;
   %ORIp(i) = p(1,2);
   
end

%% Get correlation coeffs at close range


id = find(Dist>0 & Dist<Dmax & ~isnan(dori) & ~isnan(dsf));
[r p] = corrcoef(Dist(id),dori(id))
r_oriVDist = r(1,2);
p_oriVDist = p(1,2);

[r p] = corrcoef(Dist(id),dsf(id))
r_sfVDist = r(1,2);
p_sfVDist = p(1,2);

id = find(Dist>0 & Dist<Dmax & ~isnan(dphase));
[r p] = corrcoef(Dist(id),dphase(id))
r_phaseVDist = r(1,2);
p_phaseVDist = p(1,2);


%%

figure,
subplot(4,1,1)
plot(Distdom,ORIr,'-ok'), ylabel('ori1 ori2 correlation')
xlim([0 Dmax])

subplot(4,1,2)
plot(Distdom,SFr,'-ok'),  ylabel('sf1 sf2 correlation')
idSig = find(SFp<.01);
hold on, plot(Distdom(idSig),SFr(idSig),'*r')
xlim([0 Dmax])

subplot(4,1,3)
plot(Distdom_phase,OrivSFr,'-ok'), ylabel('|dori| vs. |dsf| correlation')
xlabel('distance (mm)')
idSig = find(OrivSFp<.01);
hold on, plot(Distdom(idSig),OrivSFr(idSig),'*r')
xlim([0 Dmax])

subplot(4,1,4)
plot(Distdom_phase,Phaser,'-ok'), ylabel('phase1 phase2 correlation')
xlim([0 Dmax])

%% Plot scatter plot and error bars of pairwise difference

figure,
subplot(3,1,1)
%plot(Dist,dori,'.c')
scatter(Dist,dori,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
hold on, errorbar(Distdom,dorimu,dorisig,'LineWidth',2), ylabel('dori degrees')
hold on, plot([0 20],[1 1]*nanmean(dori),'b')
xlim([0 Dmax])
title(['r=' num2str(r_oriVDist) '; p=' num2str(p_oriVDist)])

subplot(3,1,2),
%scatter(Dist,dsf,'markerEdgeAlpha',.1,'Marker','.','LineWidth',5)
scatter(Dist,dsf,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
hold on, errorbar(Distdom,dsfmu,dsfsig,'LineWidth',2), ylabel('dsf octaves')
hold on, plot([0 20],[1 1]*nanmean(dsf),'b')
xlim([0 Dmax])
ylim([0 2])
title(['r=' num2str(r_sfVDist) '; p=' num2str(p_sfVDist)])

subplot(3,1,3)
%scatter(Dist,dphase,'Marker','o','markerEdgeColor',[1 1 1],'markerEdgeAlpha',0,...
%    'markerFaceColor',[0 0 0],'markerFaceAlpha',.05,'SizeData',50)
scatter(Dist,dphase,'.k')
hold on, errorbar(Distdom,dphasemu,dphasesig,'LineWidth',2), ylabel('dphase degrees')
hold on, plot([0 20],[1 1]*nanmean(dphase),'b')
xlim([0 Dmax])
title(['r=' num2str(r_phaseVDist) '; p=' num2str(p_phaseVDist)])
xlim([0 500])

xlabel('mean distance in Bin')
