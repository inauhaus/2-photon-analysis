function rst = orisfmap(c1,c2)

%Version 3 uses the function 'analyzeLOGtuning', which generates all of the
%sf parameters.

global cellS idExamp ACQinfo Analyzer

%Order in the looper

for p = 1:length(Analyzer.L.param)
    switch Analyzer.L.param{p}{1}
        case 'ori'
            oridim = p;
        case 's_freq'
            sfdim = p;
            
    end
end

[dum dimsort] = sort([oridim sfdim]);

%SF domains
sfdom = getdomain('s_freq');
sfdomI = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*2-1);
sfdomII = logspace(log10(sfdom(1)),log10(sfdom(end)),length(sfdom)*5-1); %strictly for plotting a smooth curve


%ORI domains
oridom = getdomain('ori');



figure
for i = 1:length(cellS.mukern)
    
    i
    kdum = cellS.mukern{i};
    kdum = permute(kdum,dimsort);  %force the kernel to have dimensions in this order
    
    
    %sforimat = squeeze(mean(kdum,1)); %mean over color
    
    
    sforimat = squeeze(kdum(:,:));
    
    oritc_r = squeeze(mean(sforimat,2)); %ori tuning curve for each repeat
    sftc_r = squeeze(mean(sforimat,1)); %sf tuning curve for each repeat
    
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
    
    
 
       
    rst.sfpref{c}(i) = params.pref;
    rst.sfhco{c}(i) = params.hco;
    rst.sfBW{c}(i) = params.BW;
    rst.sfBP{c}(i) = params.BP;
    rst.sfCoM{c}(i) = params.CoM;
    
    %rst.sftc{c}(i,:) = ffitIsf;
    rst.sftc{c}(i,:) = sftc;
    rst.sftcfit{c}(i,:) = ffitIsf;
    
    %pick a condition for the example time course
    sforimatdum = squeeze(kdum(2,:,:) + kdum(1,:,:)); %use one color to identify optimal ori and sf
    [idy idx] = find(sforimatdum == max(sforimatdum(:)));
    
    rst.tcourseOpt{c}(i,:) = squeeze(kdumT(c,idy,idx,:));
    rst.tcourseLowSF{c}(i,:) = squeeze(kdumT(c,idy,1,:));

    %plot tuning curves and fits
    if c == 1
        col = 'k';
        else
            col = 'r';
        end
        
%         if c == 1
%             figure(89)
%         else 
%             figure(90)
%         end
        figure(90)
        subplot(ceil(sqrt(length(cellS.mukern))),ceil(sqrt(length(cellS.mukern))),i)
        %         hold on
        %plot(log2(sfdomII),ffitIsf,col)
        %hold on
        %plot(log2(sfdom),sftc,['.' col])
        
        %plot(sftc_r,'.k'), 
       % hold on, 
        %plot(ffitMat,'o-b')
        %hold on,
        %plot(median(sftc_r,2),'o-r'), 
        


        %hold on
        %plot(oridom,oritc)

        h = fspecial('gaussian',size(rst.tcourseLowSF{c}(i,:)),1) - fspecial('gaussian',size(rst.tcourseLowSF{c}(i,:)),20);
        tcfilt = ifft(fft(rst.tcourseLowSF{c}(i,:)).*abs(fft(h)));
         plot(tcfilt)
         hold on
        
        %title(num2str(sfpref(i)))
        title([num2str(pfit) ' ' num2str(Gvaracc)])
        axis tight
        
  
    end
    
    %Get color selectivity e.g. %S or %S-M
    if ~isnan(rst.sfpref{1}(i)*rst.sfpref{2}(i))
        dum = mean(mean(kdum,2),3);
        %dum = [max(rst.sftc{1}(i,:)) max(rst.sftc{2}(i,:))];
        
        rst.pS(i) = dum(2)/(dum(1)+dum(2));
       
        %rst.pS(i) = rst.sftc{2}(i,1) / (rst.sftc{2}(i,1)+rst.sftc{1}(i,1));
        
        if isnan(rst.sfpref{1}(i)) | isnan(rst.sfpref{2}(i))
            rst.pS(i) = NaN;
        end

    else
        rst.pS(i) = NaN;
    end
    
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

%%
figure
scatter(abs(oridiff(rst.opref{1}*pi/180,pi/2))*180/pi,rst.pS,'.k'), 
xlabel('Ori difference from vertical; color 1 ')

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
