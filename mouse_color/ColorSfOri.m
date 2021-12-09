function rst = ColorSfOri(c1,c2)

%This was adapted from ColorSfreq.m to show some orientation tuning

global cellS idExamp ACQinfo Analyzer

%Order in the looper

for p = 1:length(Analyzer.L.param)
    switch Analyzer.L.param{p}{1}
        case 'theta'
            colordim = p;
        case 'ori'
            oridim = p;
        case 's_freq'
            sfdim = p;
            
    end
end

[dum dimsort] = sort([colordim oridim sfdim]);

sfdom = getdomain('sfreq');
sfdomI = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*2-1);

oridom = getdomain('ori');

cid = [c1 c2];

figure
for i = 1:length(cellS.mukern)   
    
    %kdum = abs(cellS.F1kern{i});
    kdum = cellS.mukern{i};
    kdum = (kdum);
    kdum = permute(kdum,dimsort);  %force the kernel to have dimensions in this order
    
    kdumsig = cellS.sigkern{i};
    kdumsig = permute(kdumsig,dimsort);  %force the kernel to have dimensions in this order
      
    kdumT = cellS.mukernTime{i};
    kdumT = permute(kdumT,[dimsort 4]);
    
    %sforimat = squeeze(mean(kdum,1)); %mean over color
    for c = 1:2 %loop each color
        
        sforimat = squeeze(kdum(cid(c),:,:));
        sforimatsig = squeeze(kdumsig(cid(c),:,:));
        
        [idy idx] = find(sforimat == max(sforimat(:)));
        
        %Compute tuning curves using a weighted average of the orthogonal
        %dimension
        sfdum = sum(sforimat); sfdum  = sfdum/sum(sfdum);
        oritc = sum(sforimat.*(ones(length(oridom),1)*sfdum),2);
        
%         oritc = mean(sforimat,2);
%         orisigtc = mean(sforimatsig,2);
        oritc = sforimat(:,idx);
        orisigtc = sforimatsig(:,idx);
        rst.oritc{c}(i,:) = oritc;
        rst.oritcsig{c}(i,:) = orisigtc;
        
        oridum = sum(sforimat,2); oridum = oridum/sum(oridum);
        sftc = sum(sforimat.*(oridum*ones(1,length(sfdom))),1);
        
   
        %Compute orientation preference and selectivity
        R = sum(phi(oritc).*exp(1i*oridom'*pi/180*2))/sum(phi(oritc));
        rst.omag{c}(i) = abs(R);
        rst.opref{c}(i) = angle(R)/2*180/pi;
        if rst.opref{c}(i)<0
            rst.opref{c}(i) = rst.opref{c}(i)+180;
        end
        
        %Compute direction preference and selectivity
        R = sum(phi(oritc).*exp(1i*oridom'*pi/180))/sum(phi(oritc));
        rst.dmag{c}(i) = abs(R);
        rst.dpref{c}(i) = angle(R)*180/pi;
        if rst.dpref{c}(i)<0
            rst.dpref{c}(i) = rst.dpref{c}(i)+360;
        end
        
        %Fit a Difference of Gaussian to the spatial frequency tuning curve
        %and id location of the peak
        sftcI = interp1(sfdom,sftc,sfdomI,'spline')';        
        [paramdum ffit_sf varacc_sf ffitIsf domIsf pk sfBW sflco sfhco] = DoGfit(sftcI,sfdomI');
        
        
        [dum sfid] = max(ffitIsf);
        rst.sfpref{c}(i) = domIsf(sfid);
        rst.sfhco{c}(i) = sfhco;
        rst.sfBW{c}(i) = sfBW;
        rst.sfBP{c}(i) = (max(sftc)-phi(sftc(1)))/(max(sftc)+phi(sftc(1)));
        
        %Compute center of mass of tuning curve
        sftcnorm = sftc - min(sftc);
        sftcnorm = sftcnorm/sum(sftcnorm);        
        rst.sfCoM{c}(i) = 2.^sum(sftcnorm.*log2(sfdom));
        
        %rst.sftc{c}(i,:) = ffitIsf;
        rst.sftc{c}(i,:) = sftc;
        rst.sftcfit{c}(i,:) = ffitIsf;
        
        %pick a condition for the example time course
        sforimatdum = squeeze(kdum(2,:,:)); %use one color for ori and sf
        [idy idx] = find(sforimatdum == max(sforimatdum(:)));
        rst.tcourse{c}(i,:) = squeeze(kdumT(c,idy,idx,:));
        idy2 = idy+length(oridom)/2;
        if idy2>length(oridom)
            idy2 = idy2-length(oridom);
        end
        rst.tcourse{c}(i,:) = (rst.tcourse{c}(i,:) + squeeze(kdumT(c,idy2,idx,:))')/2;
        
        %Remove cells with bad fits
        if varacc_sf < .6
            rst.sfpref{c}(i) = NaN;
            rst.sfhco{c}(i) = NaN;
            rst.sfBW{c}(i) = NaN;
            rst.sfBP{c}(i) = NaN;
            rst.sftc{c}(i,:) = NaN*rst.sftc{c}(i,:);
            rst.sftcfit{c}(i,:) = NaN*rst.sftcfit{c}(i,:);
            rst.opref{c}(i) = NaN;
            rst.omag{c}(i) = NaN;
            rst.dpref{c}(i) = NaN;
            rst.dmag{c}(i) = NaN;
            rst.sfCoM{c}(i) = NaN;
            
        end
               
        %plot tuning curves and fits

        subplot(ceil(sqrt(length(cellS.mukern))),ceil(sqrt(length(cellS.mukern))),i)
%         hold on
%         plot(log2(domIsf),ffitIsf)
%         hold on
%         plot(log2(sfdom),sftc,'.')
        
        plot(rst.tcourse{c}(i,:))
%         hold on
        
        %title(num2str(sfpref(i)))
        axis tight
        
  
    end
    
    %Get color selectivity e.g. %S or %S-M
    if ~isnan(rst.sfpref{1}(i)*rst.sfpref{2}(i))
        dum = mean(mean(kdum,2),3);
        %dum = [max(rst.sftc{1}(i,:)) max(rst.sftc{2}(i,:))];
        
        rst.pS(i) = dum(2)/(dum(1)+dum(2));
        
        if isnan(rst.sfpref{1}(i)) | isnan(rst.sfpref{2}(i))
            rst.pS(i) = NaN;
        end

    else
        rst.pS(i) = NaN;
    end
    
end
%%
fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
tdom =(0:(size(rst.tcourse{1},2)-1))*fr/1000; 
figure
for i = 1:length(idExamp)
     
    subplot(length(idExamp),2,(i-1)*2+1)
 
%     semilogx((sfdom),rst.sftc{1}(idExamp(i),:),'ok')
%     hold on
%     semilogx((domIsf),rst.sftcfit{1}(idExamp(i),:),'k')
%     hold on
%     semilogx((sfdom),rst.sftc{2}(idExamp(i),:),'or')
%     hold on
%     semilogx((domIsf),rst.sftcfit{2}(idExamp(i),:),'r')    
%     set(gca,'Xtick',sfdom)    

tcdum = rst.oritc{2}(idExamp(i),:)/2+rst.oritc{1}(idExamp(i),:)/2;
tcdumsig = rst.oritcsig{2}(idExamp(i),:)/2+rst.oritcsig{1}(idExamp(i),:)/2;
    errorbar(oridom,tcdum,tcdumsig,'k')    
    set(gca,'Xtick',oridom)    
    
    axis tight
    yl = ylim;
    yl(1) = 0;
    ylim(yl)
%     
ylabel('dF/F')
    subplot(length(idExamp),2,(i-1)*2+2)
    plot(tdom,rst.tcourse{1}(idExamp(i),:),'k')
    hold on
    plot(tdom,rst.tcourse{2}(idExamp(i),:),'k')
    hold on
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
    xlabel('sec'), 
    axis tight
    %title(num2str(sfpref(i)))

    
end

%%
%Normalize, then compute mean/std at each sfreq
N = length(find(~isnan(rst.sftc{1}(:,1))));
sfnorm{1} = rst.sftc{1}./(max(rst.sftc{1},[],2)*ones(1,length(sfdom)));
sfnorm{2} = rst.sftc{2}./(max(rst.sftc{2},[],2)*ones(1,length(sfdom)));
muC1 = nanmean(sfnorm{1});
muC2 = nanmean(sfnorm{2});
sigC1 = nanstd(sfnorm{1});
sigC2 = nanstd(sfnorm{2});


%% Spatial frequency comparison

figure,
subplot(4,2,1)
scatter(rst.sfpref{1},rst.sfpref{2},'.k'), hold on, plot([0 .3],[0 .3])
xlabel('Sf peak location; color 1 ')
ylabel('Sf peak location; color 2 ')
axis square
xlim([0 .3]), ylim([0 .3])

subplot(4,2,2)
histogram(log2(rst.sfpref{1}./rst.sfpref{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfpref{1}./rst.sfpref{2}));
title(['mean=' num2str(2.^mu)])

subplot(4,2,3)
scatter(rst.sfCoM{1},rst.sfCoM{2},'.k'), hold on, plot([0 .15],[0 .15])
xlabel('Sf Center of Mass; color 1 ')
ylabel('Sf Center of Mass; color 2 ')
axis square
xlim([0 .15]), ylim([0 .15])

subplot(4,2,4)
histogram(log2(rst.sfCoM{1}./rst.sfCoM{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfCoM{1}./rst.sfCoM{2}));
title(['geomean=' num2str(2.^mu)])

subplot(4,2,5)
scatter(rst.sfhco{1},rst.sfhco{2},'.k'), hold on, plot([0 .3],[0 .3])
xlabel('Sf high pass cutoff; color 1 ')
ylabel('Sf high pass cutoff; color 2 ')
axis square

subplot(4,2,6)
histogram(log2(rst.sfhco{1}./rst.sfhco{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfhco{1}./rst.sfhco{2}));
title(['mean=' num2str(2.^mu)])

subplot(4,2,7)
scatter(rst.sfBP{1},rst.sfBP{2},'.k'), hold on, plot([0 1],[0 1])
xlabel('Sf bandpass; color 1 ')
ylabel('Sf bandpass; color 2 ')
axis square


subplot(4,2,8)
histogram(rst.sfBP{1} - rst.sfBP{2},[-1:.1:1],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(rst.sfBP{1} - rst.sfBP{2});
title(['mean=' num2str(2.^mu)])

% subplot(4,2,7)
% scatter(rst.sfBW{1},rst.sfBW{2},'.k'), hold on, plot([0 4],[0 4])
% xlabel('Sf bandwidth; color 1 ')
% ylabel('Sf bandwidth; color 2 ')
% axis square


% subplot(4,2,8)
% histogram(log2(rst.sfBW{1}./rst.sfBW{2}),[-4:.5:4],'FaceColor',[1 1 1])
% set(gca,'TickDir','out')
% mu = nanmean(log2(rst.sfBW{1}./rst.sfBW{2}));
% title(['geomean=' num2str(2.^mu)])
%%

figure,errorbar(log2(sfdom),muC1/max(muC1),sigC1/sqrt(N)/max(muC1))
hold on,errorbar(log2(sfdom),muC2/max(muC2),sigC2/sqrt(N)/max(muC2))
set(gca,'XTick',log2(sfdom))
set(gca,'XTickLabel',sfdom)
xlabel('Spatial freq')
legend('color 1', 'color 2')


figure,
subplot(2,2,1), scatter(rst.pS,rst.sfpref{1}), ylabel('sf peak (color 1)'), xlabel('%color2')
subplot(2,2,2), scatter(rst.pS,rst.sfpref{2}),ylabel('sf peak (color 2)'), xlabel('%color2')
subplot(2,2,3), scatter(rst.pS,rst.sfCoM{1}), ylabel('sf CoM (color 1)'), xlabel('%color2')
subplot(2,2,4), scatter(rst.pS,rst.sfCoM{2}),ylabel('sf CoM (color 2)'), xlabel('%color2')


%% orientation comparison

figure,
subplot(2,2,1)
scatter(rst.opref{1},rst.opref{2},'.k'), hold on, plot([0 180],[0 180])
xlabel('Ori preference; color 1 ')
ylabel('Ori preference; color 2 ')
axis square

subplot(2,2,2)
scatter(rst.omag{1},rst.omag{2},'.k'), hold on, plot([0 1],[0 1])
xlabel('Ori selectivity; color 1 ')
ylabel('Ori selectivity; color 2 ')
axis square

subplot(2,2,4)
SelDiff = (rst.omag{1}-rst.omag{2})./(rst.omag{1}+rst.omag{2});
SelDiff = log2(rst.omag{1}./rst.omag{2});
histogram(SelDiff,[-3:.5:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(SelDiff);
title(['geomean=' num2str(2.^mu)])

%% direction comparison

figure,
subplot(2,2,1)
scatter(rst.dpref{1},rst.dpref{2},'.k'), hold on, plot([0 360],[0 360])
xlabel('Dir preference; color 1 ')
ylabel('Dir preference; color 2 ')
axis square

subplot(2,2,2)
scatter(rst.dmag{1},rst.dmag{2},'.k'), hold on, plot([0 1],[0 1])
xlabel('Dir selectivity; color 1 ')
ylabel('Dir selectivity; color 2 ')
axis square

subplot(2,2,4)
SelDiff = (rst.dmag{1}-rst.dmag{2})./(rst.dmag{1}+rst.dmag{2});
SelDiff = log2(rst.dmag{1}./rst.dmag{2});
histogram(SelDiff,[-3:.5:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(SelDiff);
title(['geomean=' num2str(2.^mu)])
%% Plot sf map
%%%%%%%%%%%%%%%
global maskS ACQinfo

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

cdom = jet;

% vrange = ceil(max(abs(param(:,1)))); vrange = [-vrange vrange];
% vrange = [min(param(:,1)) max(param(:,1))]

sfdomL = log2(sfdom);
sfprefL = log2(rst.sfpref{1});

vrange = [sfdomL(1) sfdomL(end)];

vertplot = (sfprefL-vrange(1))/(vrange(2)-vrange(1)) ;

vertplot = ceil(vertplot*64);
vertplot(find(vertplot>64 | vertplot<1)) = NaN;
figure
for p = 1:length(celldom) %first element is neuropil
    
    if ~isnan(vertplot(p))
        
        vc = cdom(vertplot(p),:);
        
        ax1 = subplot(2,2,1);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*xmicperpix,'.','Color',vc,'MarkerSize',30)
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
    end
    
end
colormap(ax1,cdom)

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('sfpref')
colorbar(ax1,'Ticks',[0 .5 1],'TickLabels',[sfdomL(1) 0 sfdomL(end)])
ylabel('microns')
%
% Fit a plane

%% Plot ori map
%%%%%%%%%%%%%%%

cdom = hsv;

% vrange = ceil(max(abs(param(:,1)))); vrange = [-vrange vrange];
% vrange = [min(param(:,1)) max(param(:,1))]

opref = (rst.opref{1});

vrange = [0 180];

vertplot = (opref-vrange(1))/(vrange(2)-vrange(1)) ;

vertplot = ceil(vertplot*64);
vertplot(find(vertplot>64 | vertplot<1)) = NaN;

for p = 1:length(celldom) %first element is neuropil
    
    if ~isnan(vertplot(p))
        
        vc = cdom(vertplot(p),:);
        
        ax1 = subplot(2,2,2);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*xmicperpix,'.','Color',vc,'MarkerSize',30)
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
    end
    
end
colormap(ax1,cdom)

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('opref')
colorbar(ax1,'Ticks',[0 1],'TickLabels',[0 180])
ylabel('microns')
%% Plot color map

rst.pS(find(rst.pS>1 | rst.pS<0)) = NaN;

pSmap = rst.pS(2:end);
pSmap = ceil(pSmap*64);
cdom = jet;

for p = 1:length(pSmap)
    
    if ~isnan(pSmap(p))
        
        vc = cdom(pSmap(p),:);
        
        ax2 = subplot(2,2,3);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*ymicperpix,'.','Color',vc,'MarkerSize',30)
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
    end
    
end

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('%Color2 (peak)')
vR = [0 1];
colorbar(ax2,'Ticks',[0 .5 1],'TickLabels',[0 .5 1])
colormap(ax2,cdom)

%%  Plot raw image
CHs = GetTrialData([1 0],1);
im = mean(CHs{1},3);
[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

mi = prctile(im(:),.1);
ma = prctile(im(:),99.8);
figure, subplot(1,2,1)
imagesc(xdom,ydom,im,[mi ma]), colormap gray
hold on
for i = 1:length(idExamp)    
   plot(CoM(idExamp(i)-1,2)*xmicperpix,CoM(idExamp(i)-1,1)*ymicperpix,'or','MarkerSize',12)     
end
axis image
subplot(1,2,2)
imagesc(xdom,ydom,maskS.im{1}), colormap gray
hold on
for i = 1:length(idExamp)    
   plot(CoM(idExamp(i)-1,2)*xmicperpix,CoM(idExamp(i)-1,1)*ymicperpix,'or','MarkerSize',12)     
end
axis image
title('local x-correlation')
