function BPFvColor(Q,idyield,SOA_examp,SOB_examp,DO_examp,Lum_examp)


idB = idyield;

cpref = (Q.Colpk-Q.Lumpk)./(Q.Colpk+Q.Lumpk);

x = cpref(idB);
y = Q.tc.x45_135.sfBPB(idB);
figure
subplot(4,1,1)
scatter(x,y,'.k'), hold on 
scatter(x(DO_examp),y(DO_examp),'or')
scatter(x(SOA_examp),y(SOA_examp),'ob')
scatter(x(SOB_examp),y(SOB_examp),'ok')
scatter(x(Lum_examp),y(Lum_examp),'og')
hold on, plot([-1 1],[.5 .5],'--k'), plot([0 0],[0 1],'--k')
ylim([-.05 1.05]), xlim([-1.05 1.05]), 
xlabel('%color'),ylabel('sf BPF (S-R)'),
[r p] = corrcoef(x,y)
r = r(1,2); p = p(1,2);
title(['r=' num2str(r) '; p=' num2str(p)])
DOid = find(x>0 & y>0.5)
N_doubleOpponent = length(DOid)
getBPFvColQuadN(x, y, 1)
axis square

x = cpref(idB);
y = Q.tc.x45_135.sfBPA(idB);
subplot(4,1,2)
scatter(x,y,'.k'), hold on 
scatter(x(DO_examp),y(DO_examp),'or')
scatter(x(SOA_examp),y(SOA_examp),'ob')
scatter(x(SOB_examp),y(SOB_examp),'ok')
scatter(x(Lum_examp),y(Lum_examp),'og')
hold on, plot([-1 1],[.5 .5],'--k'), plot([0 0],[0 1],'--k')
ylim([-.05 1.05]), xlim([-1.05 1.05]), 
xlabel('%color'),ylabel('sf BPF (S+R)'),
[r p] = corrcoef(x,y)
r = r(1,2); p = p(1,2);
title(['r=' num2str(r) '; p=' num2str(p)])
DOid = find(x>0 & y>0.5)
N_doubleOpponent = length(DOid)
getBPFvColQuadN(x, y, 1)
axis square

x = cpref(idB);
y = Q.tc.x0_90.sfBPB(idB);
subplot(4,1,3)
scatter(x,y,'.k'), hold on 
scatter(x(DO_examp),y(DO_examp),'or')
scatter(x(SOA_examp),y(SOA_examp),'ob')
scatter(x(SOB_examp),y(SOB_examp),'ok')
scatter(x(Lum_examp),y(Lum_examp),'og')
hold on, plot([-1 1],[.5 .5],'--k'), plot([0 0],[0 1],'--k')
ylim([-.05 1.05]), xlim([-1.05 1.05]), 
xlabel('%color'), ylabel('sf BPF (S)')
[r p] = corrcoef(x,y)
r = r(1,2); p = p(1,2);
title(['r=' num2str(r) '; p=' num2str(p)])
DOid = find(x>0 & y>0.5)
N_doubleOpponent = length(DOid)
getBPFvColQuadN(x, y, 1)
axis square

x = cpref(idB);
y = Q.tc.x0_90.sfBPA(idB);
subplot(4,1,4)
scatter(x,y,'.k'), hold on 
scatter(x(DO_examp),y(DO_examp),'or')
scatter(x(SOA_examp),y(SOA_examp),'ob')
scatter(x(SOB_examp),y(SOB_examp),'ok')
scatter(x(Lum_examp),y(Lum_examp),'og')
hold on, plot([-1 1],[.5 .5],'--k'), plot([0 0],[0 1],'--k')
ylim([-.05 1.05]), xlim([-1.05 1.05]), 
xlabel('%color'),,ylabel('sf BPF (R)')
[r p] = corrcoef(x,y)
r = r(1,2); p = p(1,2);
title(['r=' num2str(r) '; p=' num2str(p)])
DOid = find(x>0 & y>0.5)
N_doubleOpponent = length(DOid)
getBPFvColQuadN(x, y, 1)
axis square
