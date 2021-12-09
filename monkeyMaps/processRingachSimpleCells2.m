load('C:\2pScanboxAnalysis\RF_db.mat')
pF0


%%
figure,
clear Allparam Allffit Allvaraccount AllparamGuess
k = 1;
for i = 1:10
    for j = 1:25
        
         G = getGaborGuess(rf(k).rf,rf(k).d2p);
         AllparamGuess(k,:) = [G.ymu G.xmu G.ysig G.xsig G.ori G.sf G.phase];
%         
%         %Convert units to degrees
         pixperim = size(rf(k).rf,1);
         AllparamGuess(k,1:4) = AllparamGuess(k,1:4)/rf(k).d2p;
         AllparamGuess(k,6) = AllparamGuess(k,6)*rf(k).d2p/pixperim; %Convert cycles/image to cyc/degree

        [Allparam(k,:) Allffit{k} Allvaraccount(k)] = Gaborfit2Drot(rf(k).rf,rf(k).d2p);
        
   
        subplot(10,25,k)
        
        imagesc(Allffit{k})
        axis off
        
        nx = 4*Allparam(k,4)*Allparam(k,6);
        title(num2str(round(nx*1000)/1000))
        
        k = k+1;
        drawnow
    end
end

%% Fit Gaussian in ori domain

sfmodel  = 'LogGauss';
%sfmodel  = 'LinGauss';

figure,
%clear Allparam Allffit Allvaraccount AllparamGuess
k = 1;
for i = 1:10
    for j = 1:25
        
%         
%         %Convert units to degrees
         pixperim = size(rf(k).rf,1);
%          AllparamGuess(k,1:4) = AllparamGuess(k,1:4)/rf(k).d2p;
%          AllparamGuess(k,6) = AllparamGuess(k,6)/pixperim*rf(k).d2p; %Convert cycles/image to cyc/degree

        rfF = fftshift(fftshift(abs(fft2(rf(k).rf)),1),2);
        rfF = rfF/max(rfF(:));
        
       
        %[AllparamF(k,:) AllffitF{k} AllvaraccountF(k)] = Gaussfit2Dpolar(rfF);
%         [param ffit varacc sigma] = oriFitfrom2D(rfF);
%         oriSig1D(k) = sigma;        
%         varacc1Dori(k) = varacc;
%         
%         
%         [param  varacc pref sig sfMarg domMarg ffit domII] = sfFitfrom2D(rfF,rf(k).d2p);
        
        [sfpref sfsig sfsigOct orisig varacc_sf varacc_ori] = sforiFitfrom2D_2(rfF,rf(k).d2p,0,sfmodel);
        
        
        oriSig1D_raw(k) = orisig;
        sfSig1D_raw(k) = sfsig; 
        sfSigOct1D_raw(k) = sfsigOct;
        sfpref1D_raw(k) = sfpref;
        varacc1Dsf_raw(k) = varacc_sf;
        varacc1Dori_raw(k) = varacc_ori;
        
        
        
        rfF = fftshift(fftshift(abs(fft2(Allffit{k})),1),2);
        rfF = rfF/max(rfF(:));
        
%          [sfpref sfsig sfsigOct orisig varacc_sf varacc_ori] = sforiFitfrom2D_2(rfF,rf(k).d2p,0,sfmodel);
%         
%         
%         oriSig1D_fit(k) = orisig;
%         sfSig1D_fit(k) = sfsig; 
%         sfSigOct1D_fit(k) = sfsigOct; 
%         sfpref1D_fit(k) = sfpref;
%         varacc1Dsf_fit(k) = varacc_sf;
%         varacc1Dori_fit(k) = varacc_ori; 
        
  
        subplot(10,25,k) 
        
%         ymarg = mean(rfF,2);
%         id = find(ymarg<.0001);
%         rfF(id,:) = [];
%         
%         xmarg = mean(rfF,1);
%         id = find(xmarg<.0001);
%         rfF(:,id) = [];

        size(rfF)
        
        %plot(domMarg,sfMarg), hold on, plot(domII,ffit)
        %xlim([domMarg(1) domMarg(end)])
        imagesc(rfF)
        axis off
        axis square
        
        %orisig = AllparamF(k,2);
        title(num2str(round(orisig)))
        %title(num2str(round(sfpref1D(k)*10)/10))
        %title(num2str(varacc1Dori(k)))
        
        k = k+1;
        drawnow
    end
end
%%
sfpref_Gabor = Allparam(:,6);

figure,
subplot(2,2,1)
id = find(Allvaraccount(:)>.4 & Allparam(:,6)>.1 & Allparam(:,4)>.01);
loglog(sfpref_Gabor(id),oriSig1D_fit(id),'.k')
xlabel('Gabor fit SF'), ylabel('ori sig of fft(fit)')
[r p] = corrcoef(log2(sfpref_Gabor(id)),log2(oriSig1D_fit(id)))
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])
subplot(2,2,3)
loglog(sfpref_Gabor(id),oriSig1D_raw(id),'.k')
xlabel('Gabor fit SF'), ylabel('ori sig of fft(raw)')
[r p] = corrcoef(log2(sfpref_Gabor(id)),log2(oriSig1D_raw(id)))
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])
subplot(2,2,2)
loglog(sfpref1D_fit(id),oriSig1D_fit(id),'.k')
xlabel('SF pref of fft(fit)'), ylabel('ori sig of fft(fit)')
[r p] = corrcoef(log2(sfpref1D_fit(id)),log2(oriSig1D_fit(id)))
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])
subplot(2,2,4)
loglog(sfpref1D_fit(id),oriSig1D_raw(id),'.k')
[r p] = corrcoef(log2(sfpref1D_fit(id)),log2(oriSig1D_raw(id)))
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])
xlabel('SF pref of fft(fit)'), ylabel('ori sig of fft(raw)')


%%

id = find(Allvaraccount(:)>.4 & Allparam(:,6)>.1 & Allparam(:,4)>.01);
%id = 1:size(AllparamGuess,1);
%param = AllparamGuess;

%id = find(Allvaraccount(:)>.5 & Allparam(:,4)>.01 & sfpref1D'>.1 & varacc1Dsf' > .6);

param = Allparam;

%ypos = param(id,1);
%xpos = param(id,2);
ysig = param(id,3);
xsig = param(id,4);
sfpref = param(id,6);
%sfpref = AllparamGuess(id,6);
%sfpref = sfpref1D(id)';
oripref = param(id,3);
%orisig = oriSig1D_raw(id);

a = (1./(4*sfpref));
b = (xsig);

%a = log2(1./(4*sfpref));
%b = log2(xsig);
[param ffit varaccount dom] = NLsoftsat(a,b,[1 .5]);
%%
%b = ysig./xsig;
%b = log2(1./(ysig.*sfpref*4));
figure,loglog(a,b,'.k'),
[r p] = corrcoef(a,b);

%sfpref = AllparamGuess(id,6);
% a2p = (1./(4*TCroAll.sfpref));
% b2p = (TCopAll.profileSize/2); %actual RF width; sigma
% %b2p = TCopAll.ysize*2;
% 
% id = find(isnan(a2p.*b2p));
% b2p(id) = [];
% a2p(id) = [];

%a2p = log2(a2p);
%b2p = log2(b2p);


xlabel('0.25/(SF preference)'), ylabel('RF size; sigma'),% xlim([0 3])
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])
axis square
%[param ffit varaccount] = Sigfit(a,b)



%alpha in paper is param(1)/4
%T in paper is param(2)

%[param ffit dom varaccount] = Saturationfit2((a),(b),[1 .1 1])

%[param ffit varaccount domu] = Expfit2(a,b,[2 1 0])
% yint = param(3)+param(1)
%[param ffit varaccount dom] = Decay_fit(b,a,[.2 0 1]);
%yint = 1/param(1) + param(2)
% hold on
% plot(dom,ffit,'r')
%set(gca,'XTick',[.5 1 2 4 8],'YTick',[.25 .5 1 2])


hold on
loglog(dom,dom,'k')

hold on, loglog(dom,ffit,'--k')

set(gca,'XTick',[1/16 1/8 .25 .5 1 2]/2,'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
set(gca,'YTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlim([1/32 1])
ylim([1/32 1])


% hold on,
% loglog(a2p,b2p,'.k')

legend('Ringach data','Ringach fit','4sig = 1/sf','2p data')
%axis square

%% Aspect ratio (ori selectivity) vs. spatial frequency

%id = find(Allvaraccount > .4 &  varacc1Dori > .7 & Allparam(:,6)'>.1);

%id = find(varacc1Dsf > .7 &  varacc1Dori > .7);

%Ringach simple cells
%param = Allparam;
%sfpref = Allparam(id,6); %carrier fit in spatial domain
sfpref = sfpref1D(id);  %1D fit from radial marginal
sfpref = AllparamGuess(id,6); %peak of smoothed sfdomin
orisig = oriSig1D(id);
a = log2(sfpref);
b = orisig;

%2p data
b2p = TCroAll.orisig;
a2p = log2(TCroAll.sfpref);
id = find(isnan(a2p.*b2p));
b2p(id) = [];
a2p(id) = [];


%b = atan(ysig.*sfpref)*180/pi; %ori sigma from spatial domain

%b = ysig./xsig;
%b = log2(1./(ysig.*sfpref*4));
figure,
subplot(2,2,1)
plot(a,b,'ok'),
[r1 p1] = corrcoef((a),b);
[r2 p2] = corrcoef((a2p),b2p);
xlim([-2 3])

%ylim([0 2.5]), xlim([0 7])

xlabel('spatial frequency preference'), ylabel('orientation bandwidth (sig)'),% xlim([0 3])
title(['r(Dario)=' num2str(r1(1,2)) '; p=' num2str(p1(1,2)) '; r(2p)=' num2str(r2(1,2)) '; p=' num2str(p2(1,2))])

%set(gca,'XTick',[.5 1 2 4 8])
%set(gca,'YTick',[.25 .5 1 2])
ylim([0 90])


hold on,
plot(a2p,b2p,'ob')

legend('Ringach data','2p data')
axis square


%histogram of ori bandwidths
bins = 0:2:100; 
[h edges] = histcounts(b,bins); 
subplot(2,2,2)
plot(edges(2:end),cumsum(h/sum(h)),'k')
hold on
[h edges] = histcounts(b2p,bins),
plot(edges(2:end),cumsum(h/sum(h)),'b')
xlim([0 45])


legend('Ringach','2p')
xlabel('ori bandwidth (sigma)')
title(['Ringach med = ' num2str(median(b)) '; 2p med = ' num2str(median(b2p)) ])

%histogram of sf prefs
%bins = logspace(log10(.1),log10(8),10);
bins = linspace(log2(.1),log2(8),10);
[h edges] = histcounts(a,bins); 
subplot(2,2,3)
plot(edges(2:end),cumsum(h/sum(h)),'k')
hold on
[h edges] = histcounts(a2p,bins),
plot(edges(2:end),cumsum(h/sum(h)),'b')
%xlim([0 45])
legend('Ringach','2p')
xlabel('SF preference (octaves)')
title(['Ringach geomean = ' num2str(2.^mean(a)) '; 2p med = ' num2str(2.^mean(a2p)) ])
%xlim([bins(2) bins(end)])
xlim([-2 3])
%%
a = log2((1./sfpref));
b = log2((4*xsig));

figure,scatter(a,b)
hold on, plot([-2 2],[-2 2]),
xlabel('1/sf'), 
ylabel('4sig')
axis square 

PercBins = [0:20:100];
clear muPeriod muWidth sigWidth
for i = 1:length(PercBins)-1
   
    BinLow = prctile(a,PercBins(i));
    BinHigh = prctile(a,PercBins(i+1));
    idBin = find(a > BinLow  &  a <= BinHigh);
    muPeriod(i) = mean(a(idBin));
    muWidth(i) = median(b(idBin));
    sigWidth(i) = std(b(idBin));
    
end

[param ffit dom varaccount] = Saturationfit2(a,b,[1 1 0])

hold on, plot(dom,ffit)

hold on
errorbar(muPeriod,muWidth,sigWidth,'k')

 
 %%
 k = 1;
 figure
radonDom = 0:1:179;
hh = fspecial('gaussian',size(magim),5);
 for i = 1:100
     RF = rf(i).rf;
     magim = RF.^2;
     hh = fspecial('gaussian',size(magim),4);
     LocRMS = sqrt(ifft2(fft2(magim).*abs(fft2(hh))));
     
     [idy idx] = find(LocRMS == max(LocRMS(:)));
     dim = size(LocRMS);
     
     LocRMS = circshift(LocRMS,[dim(1)/2-idy dim(2)/2-idx]);
     RF = circshift(RF,[dim(1)/2-idy dim(2)/2-idx]);
     RF = RF.*(hann(dim(1))*hann(dim(2))');
 
     subplot(10,10,k)
     k = k+1;
    imagesc(RF)
    axis image
    title('hi')
    axis off
     
     
 end
     

