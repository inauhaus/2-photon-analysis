function rst = ColorSfreq4(c1,c2)

%Version 3 uses the function 'analyzeLOGtuning', which generates all of the
%sf parameters.

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

%SF domains
sfdom = getdomain('s_freq');
sfdomI = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*2-1);
sfdomII = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*5-1); %strictly for plotting a smooth curve

%These are used for plotting outside of this function
rst.sfdom = sfdom;
rst.sfdomI = sfdomI;
rst.sfdomII = sfdomII;

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

figure
for i = 1:length(cellS.mukern)   

   i
    kdum = cellS.mukern{i};
    kdum = permute(kdum,dimsort);  %force the kernel to have dimensions in this order
    
    kdum_r = squeeze(cellS.Repkern{i});
    kdum_r = permute(kdum_r,[dimsort 4]);  %force the kernel to have dimensions in this order.  Last dimensions is repeats    
      
    kdumT = squeeze(cellS.mukernTime{i});
    kdumT = permute(kdumT,[dimsort 4]);
    
    %sforimat = squeeze(mean(kdum,1)); %mean over color
    for c = 1:2 %loop each color
        
        sforimat = squeeze(kdum(cid(c),:,:));
        rst.sforimat{c}(:,:,i) = squeeze(mean(sforimat,3));
        
        sforimat_r = squeeze(kdum_r(cid(c),:,:,:));
        oritc_r = squeeze(mean(sforimat_r,2)); %ori tuning curve for each repeat
        sftc_r = squeeze(mean(sforimat_r,1)); %sf tuning curve for each repeat
        
        oritc = mean(oritc_r,2); %average over repeats
        sftc = mean(sftc_r,2);
        
        %Compute orientation preference and selectivity non-parametrically
        R = sum(phi(oritc).*exp(1i*oridom'*pi/180*2))/sum(phi(oritc));
        rst.omag{c}(i) = abs(R);  %non-parametric selectivity
        rst.opref{c}(i) = angle(R)/2*180/pi;
        if rst.opref{c}(i)<0
            rst.opref{c}(i) = rst.opref{c}(i)+180;
        end
        rst.oritc{c}(i,:) = oritc;
        
        %Fit a Gaussian to ori tuning curve to get bandwidth
        %oritcI = interp1(oridom,oritc,oridomI,'spline')';        
        %[paramdum ffit_sf varacc_sf ffitIsf domIsf pk sfBW sflco sfhco] = DoGfit(sftcI,sfdomI');
        
        %[mu sigma ffit varacc Gparam] = CircGaussFit(oritcI');
        %Gparam(2) = sigma;
        %[Gparam ffit Gvaracc] = Gaussfit(oridomI,oritcI',1);
        
        %rst.osig{c}(i)  = Gparam(2);
                
        %ffitIori = interp1(oridomI,ffit,oridomII,'spline')';  %This is just used for plotting
        
        %rst.oritcfit{c}(i,:) = ffitIori;
        
        %Compute direction preference and selectivity
        R = sum(phi(oritc).*exp(1i*oridom'*pi/180))/sum(phi(oritc));
        rst.dmag{c}(i) = abs(R);
        rst.dpref{c}(i) = angle(R)*180/pi;
        if rst.dpref{c}(i)<0
            rst.dpref{c}(i) = rst.dpref{c}(i)+360;
        end
        
        %Fit a Gaussian to the spatial frequency tuning curve
        %and id location of the peak

        cothresh = 0.75; % high pass cutoff threshold: percent of max peak

        %[params ffitIsf domI domII pfit Gvaracc repFitQuality] = analyzeLOGtuning(sftc(:)',sfdom,cothresh,sftc_r);
        
        [params ffitIsf domI domII pfit Gvaracc repFitQuality] = analyzeLINtuning_DoG2(sftc(:)',sfdom,cothresh,sftc_r);
        
        [params.lco params.hco] = gethcolco(domII,ffitIsf',cothresh); 
        params.BW = params.hco - params.lco; %2sigma bandwidth in cyc/deg from Gaussian fit  
        
        %Regenerate the fit w/o downsampling
        ffit = exp(-(log2(sfdom)-params.raw(1)).^2/(2*params.raw(2).^2));
        %ffit = params.raw(3)*ffit + params.raw(4);
        ffit = params.raw(3)*ffit ;
        ffitMat = ffit(:)*ones(1,getnorepeats(1));
        varacc = (var(sftc_r(:))-var(sftc_r(:)-ffitMat(:)))/var(sftc_r(:));
        [dum pfit] = corrcoef(sftc_r(:),ffitMat(:));
        pfit = pfit(1,2);
 
        rst.sfpref{c}(i) = params.pref;       
        rst.sfhco{c}(i) = params.hco;
        rst.sfBW{c}(i) = params.BW;
        rst.sfBP{c}(i) = params.BP;
        rst.sfCoM{c}(i) = params.CoM;       
        
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
        if Gvaracc < .7
        %if pfit > .05
        %if repFitQuality.p > .05
            rst.sfpref{c}(i) = NaN;
            rst.sfhco{c}(i) = NaN;
            rst.sfBW{c}(i) = NaN;
            rst.sfBP{c}(i) = NaN;
            rst.sftc{c}(i,:) = NaN*rst.sftc{c}(i,:);
            rst.sftcfit{c}(i,:) = NaN*rst.sftcfit{c}(i,:);
            rst.oritc{c}(i,:) = NaN*rst.oritc{c}(i,:);
            rst.opref{c}(i) = NaN;
            rst.omag{c}(i) = NaN;
            rst.osig{c}(i) = NaN;
            rst.dpref{c}(i) = NaN;
            rst.dmag{c}(i) = NaN;
            rst.sfCoM{c}(i) = NaN;            
            
        end
               
        %plot tuning curves and fits
        if c == 1
            col = 'k';
        else
            col = 'r';
        end
        
        if c == 1
            figure(89)
        else 
            figure(90)
        end
        subplot(ceil(sqrt(length(cellS.mukern))),ceil(sqrt(length(cellS.mukern))),i)
        %         hold on
        %plot(log2(sfdomII),ffitIsf,col)
        %hold on
        %plot(log2(sfdom),sftc,['.' col])
        
        plot(sftc_r,'.k'), 
       % hold on, 
        %plot(ffitMat,'o-b')
        hold on,
        plot(median(sftc_r,2),'o-r'), 
        


        %hold on
        %plot(oridom,oritc)

        
%         plot(rst.tcourse{c}(i,:))
%         hold on
        
        %title(num2str(sfpref(i)))
        title([num2str(pfit) ' ' num2str(Gvaracc)])
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
oridomplot = [oridom 180];


fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
tdom =(0:(size(rst.tcourse{1},2)-1))*fr/1000; 
figure
for i = 1:length(idExamp)
    
    subplot(length(idExamp),5,(i-1)*5+1)
    imagesc(sfdom,oridom,rst.sforimat{1}(:,:,idExamp(i))), axis square
    set(gca,'YTick',[0 180])
    if i == length(idExamp)
        xlabel('SF (cyc/deg)')
        ylabel('Direction (deg)')
    end
    if i == 1
       title('color 1') 
    end
    
    subplot(length(idExamp),5,(i-1)*5+2)
    imagesc(sfdom,oridom,rst.sforimat{2}(:,:,idExamp(i))), axis square
    set(gca,'YTick',[0 180])
    colormap gray    
    if i == 1
       title('color 2') 
    end
     
    subplot(length(idExamp),5,(i-1)*5+3)
    
    polarplot([oridom 360]*pi/180,[rst.oritc{1}(idExamp(i),:) rst.oritc{1}(idExamp(i),1)],'.-k')
    hold on
    polarplot([oridom 360]*pi/180,[rst.oritc{2}(idExamp(i),:) rst.oritc{2}(idExamp(i),1)],'.-r')
    hold on,
    polarplot([0 0],[0 0],'.k')
    axis off
    %ax = gca;
    %ax.ThetaTick = [0 90 180 270];
    %set(gca,'ThetaTick',[0 90 180 270])
    %hold on
    %plot([oridomII 180],[rst.oritcfit{2}(idExamp(i),:) rst.oritcfit{2}(idExamp(i),1)],'r')  
%     xlim([0 360])
%     axis tight
%     yl = ylim;
%     yl(1) = min([0 yl(1)]);
%     ylim(yl)
%     set(gca,'XTick',[0 90 180])
    
    
    subplot(length(idExamp),5,(i-1)*5+4)
 
    semilogx((sfdom),rst.sftc{1}(idExamp(i),:),'ok')
    hold on
    semilogx((domII),rst.sftcfit{1}(idExamp(i),:),'k')
    hold on
    semilogx((sfdom),rst.sftc{2}(idExamp(i),:),'or')
    hold on
    semilogx((domII),rst.sftcfit{2}(idExamp(i),:),'r')    
    set(gca,'Xtick',sfdom)    
    
    axis tight
    yl = ylim;
    yl(1) = min([0 yl(1)]);
    ylim(yl)
%     
    subplot(length(idExamp),5,(i-1)*5+5)
    plot(tdom,rst.tcourse{1}(idExamp(i),:),'k')
    hold on
    plot(tdom,rst.tcourse{2}(idExamp(i),:),'r')
    hold on
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
    xlabel('sec'), ylabel('dF/F')
    axis tight
    %title(num2str(sfpref(i)))

    
end


%% Spatial frequency comparison

figure,
subplot(4,2,1)
loglog(rst.sfpref{1},rst.sfpref{2},'.k'), hold on, loglog([sfdom(1) .3],[sfdom(1) .3],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf peak loc; color 1 ')
ylabel('Sf peak loc; color 2 ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])

subplot(4,2,2)
histogram(log2(rst.sfpref{1}./rst.sfpref{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfpref{1}./rst.sfpref{2}));
[h p] = ttest(log2(rst.sfpref{1}./rst.sfpref{2}));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,3)
loglog(rst.sfCoM{1},rst.sfCoM{2},'.k'), hold on, loglog([sfdom(1) .15],[sfdom(1) .15],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf Center of Mass; color 1 ')
ylabel('Sf Center of Mass; color 2 ')
axis square
xlim([0 .15]), ylim([0 .15])

subplot(4,2,4)
histogram(log2(rst.sfCoM{1}./rst.sfCoM{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfCoM{1}./rst.sfCoM{2}));
[h p] = ttest(log2(rst.sfCoM{1}./rst.sfCoM{2}));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,5)
loglog(rst.sfhco{1},rst.sfhco{2},'.k'), hold on, loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf high pass cutoff; color 1 ')
ylabel('Sf high pass cutoff; color 2 ')
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
axis square

subplot(4,2,6)
histogram(log2(rst.sfhco{1}./rst.sfhco{2}),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(rst.sfhco{1}./rst.sfhco{2}));
[h p] = ttest(log2(rst.sfhco{1}./rst.sfhco{2}));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,7)
scatter(rst.sfBP{1},rst.sfBP{2},'.k'), hold on, plot([0 1],[0 1])
xlabel('Sf bandpass; color 1 ')
ylabel('Sf bandpass; color 2 ')
axis square

subplot(4,2,8)
histogram(rst.sfBP{1} - rst.sfBP{2},[-1:.1:1],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmedian(rst.sfBP{1} - rst.sfBP{2});
[h p] = ttest(rst.sfBP{1} - rst.sfBP{2});
title(['mu=' num2str(mu) ' p=' num2str(p)])
xlabel('x-y')


%%

%Normalize, then compute mean/std at each sfreq
N = length(find(~isnan(rst.sftc{1}(:,1))));
sfnorm{1} = rst.sftc{1}./(max(rst.sftc{1},[],2)*ones(1,length(sfdom)));
sfnorm{2} = rst.sftc{2}./(max(rst.sftc{2},[],2)*ones(1,length(sfdom)));
muC1 = nanmean(sfnorm{1});
muC2 = nanmean(sfnorm{2});
sigC1 = nanstd(sfnorm{1});
sigC2 = nanstd(sfnorm{2});

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
% subplot(2,2,1)
% scatter(rst.opref{1},rst.opref{2},'.k'), hold on, plot([0 180],[0 180])
% xlabel('Ori preference; color 1 ')
% ylabel('Ori preference; color 2 ')
% axis square

subplot(2,2,1)
scatter(abs(oridiff(rst.opref{1}*pi/180,pi/2))*180/pi,abs(oridiff(rst.opref{2}*pi/180,pi/2))*180/pi,'.k'), hold on, plot([0 90],[0 90])
xlabel('Ori difference from horizontal; color 1 ')
ylabel('Ori difference from horizontal; color 2 ')
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

%%%%%%SF map from luminance%%%%%%%%%

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
title('sfpref; Color 1')
colorbar(ax1,'Ticks',[0 1],'TickLabels',[sfdom(1) sfdom(end)])
ylabel('microns')

%%%%%%SF map from Color%%%%%%%%%%%%%%%%%%%%%%%%%

sfprefC = log2(rst.sfpref{2});

vrange = [sfdomL(1) sfdomL(end)];

vertplot = (sfprefC-vrange(1))/(vrange(2)-vrange(1)) ;

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
title('sfpref; Color 2')
colorbar(ax1,'Ticks',[0 1],'TickLabels',[sfdom(1) sfdom(end)])
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
        
        ax1 = subplot(2,2,3);
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
SMmap = jet;
SMmap(:,1) = .4;
SMmap = SMmap(1:40,:);
SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
SMmapdum = SMmap;
SMmap(:,2) = SMmapdum(:,3);
SMmap(:,3) = SMmapdum(:,2);


rst.pS(find(rst.pS>1 | rst.pS<0)) = NaN;

pSmap = rst.pS(2:end);
pSmap = ceil(pSmap*64);
cdom = SMmap;

for p = 1:length(pSmap)
    
    if ~isnan(pSmap(p))
        
        vc = cdom(pSmap(p),:);
        
        ax2 = subplot(2,2,4);
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
try
    im = 0;
    Ntrials = 5;
    for i = 1:Ntrials
        CHsdum = GetTrialData([1 0],1);
        im = im+median(CHsdum{1},3);
    end
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
catch
    
    'Could not load raw image. No access to raw file'
    
end
