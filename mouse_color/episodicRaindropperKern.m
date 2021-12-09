pF0

%initialize the Gui

global G_handles Analyzer cellS maskS

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-400');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%% Load raindropper expt
anim = 'ny6';

expt = 'u003_006'; %use mask from other expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u003_009'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 3_6'];
load(tracepath,'cellS')
getCellStats

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories




global cellS maskS


xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;

KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

'e:\2p_data\ny6\ny6_003_009'
'c:\AnalyzerFiles\ny6\ny6_u003_009'
%Celltraces called: 'ny6_u003_009_cellS_aligned to 3_6'
%Same mask as 3_6'

'e:\2p_data\ny6\ny6_002_006'
'c:\AnalyzerFiles\ny6\ny6_u002_006'

starttime = 100;
stoptime = 1200;

basestart = -200;
basestop = 0;

%%

CR = 3; %Contrast ratio

clear rst
kern = cellS.mukern;

figure
for pid = 2:length(kern)
    p = pid-1;
    
    kern{pid}(1) = kern{pid}(1)*CR;
    kern{pid}(2) = kern{pid}(2)*CR;
    
    rst.sOFF.kern(p,:) = CR*kern{pid}(4,:)';
    rst.sON.kern(p,:) = CR*kern{pid}(2,:)';
    rst.mOFF.kern(p,:) = kern{pid}(3,:)';
    rst.mON.kern(p,:) = kern{pid}(1,:)';
    
    rst.Env.kern(p,:) = mean(kern{pid},1)';
    rst.EnvON.kern(p,:) = rst.sON.kern(p,:) + rst.mON.kern(p,:);
    rst.EnvOFF.kern(p,:) = rst.sOFF.kern(p,:) + rst.mOFF.kern(p,:);
    rst.EnvM.kern(p,:) = rst.mON.kern(p,:) + rst.mOFF.kern(p,:);
    rst.EnvS.kern(p,:) = rst.sON.kern(p,:) + rst.sOFF.kern(p,:);
    
  
    [rst.sOFF.param(p,:) rst.sOFF.ffit(p,:) rst.sOFF.varacc] = Gaussfit(posdom,rst.sOFF.kern(p,:),0);
    [rst.sON.param(p,:) rst.sON.ffit(p,:) rst.sON.varacc] = Gaussfit(posdom,rst.sON.kern(p,:),0);
    [rst.mOFF.param(p,:) rst.mOFF.ffit(p,:) rst.mOFF.varacc] = Gaussfit(posdom,rst.mOFF.kern(p,:),0);
    [rst.mON.param(p,:) rst.mON.ffit(p,:) rst.mON.varacc] = Gaussfit(posdom,rst.mON.kern(p,:),0);
    
    [rst.Env.param(p,:) rst.Env.ffit(p,:) rst.Env.varacc] = Gaussfit(posdom,rst.Env.kern(p,:),0);
    [rst.EnvON.param(p,:) rst.EnvON.ffit(p,:) rst.EnvON.varacc] = Gaussfit(posdom,rst.EnvON.kern(p,:),0);
    [rst.EnvOFF.param(p,:) rst.EnvOFF.ffit(p,:) rst.EnvOFF.varacc] = Gaussfit(posdom,rst.EnvOFF.kern(p,:),0);
    [rst.EnvM.param(p,:) rst.EnvM.ffit(p,:) rst.EnvM.varacc] = Gaussfit(posdom,rst.EnvM.kern(p,:),0);
    [rst.EnvS.param(p,:) rst.EnvS.ffit(p,:) rst.EnvS.varacc] = Gaussfit(posdom,rst.EnvS.kern(p,:),0);

    
    figure(99)
    subplot(10,8,p)
    plot(posdom,rst.sOFF.ffit(p,:),'b'), hold on, plot(posdom,rst.sOFF.kern(p,:),'.b')
    hold on
    plot(posdom,rst.sON.ffit(p,:),'r'), hold on, plot(posdom,rst.sON.kern(p,:),'.r')
    if p == 1; title('SOFF SON'), end
    
    figure(100)
    subplot(10,8,p)
    plot(posdom,rst.mOFF.ffit(p,:),'b'), hold on, plot(posdom,rst.mOFF.kern(p,:),'.b')
    hold on
    plot(posdom,rst.mON.ffit(p,:),'r'), hold on, plot(posdom,rst.mON.kern(p,:),'.r')
    if p == 1; title('MOFF MON'), end
    
    figure(101)
    subplot(10,8,p)
    plot(posdom,rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:),'g'), hold on, 
    hold on
    plot(posdom,rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:),'b'),
    if p == 1; title('S & M'), end
    
    figure(102)
    subplot(10,8,p)
    plot(posdom,rst.mOFF.ffit(p,:),'g'), hold on, plot(posdom,rst.mOFF.kern(p,:),'.g')
    hold on
    plot(posdom,rst.sOFF.ffit(p,:),'b'), hold on, plot(posdom,rst.sOFF.kern(p,:),'.b')
    if p == 1; title('MOFF SOFF'), end
    
    figure(103)
    subplot(10,8,p)
    plot(posdom,rst.mON.ffit(p,:),'g'), hold on, plot(posdom,rst.mON.kern(p,:),'.g')
    hold on
    plot(posdom,rst.sON.ffit(p,:),'b'), hold on, plot(posdom,rst.sON.kern(p,:),'.b')
    if p == 1; title('MON SON'), end
    
 
    
%     Emat(p,:) = sqrt(mean(kern{pid}.^2,2));    
%     Emat2(p,1) = param(p,3)*param(p,2);
%     Emat2(p,2) = paramS(p,3)*paramS(p,2);
%     Emat2(p,3) = paramM(p,3)*paramM(p,2);
    %plot(kern{p}')
    
    if rst.sOFF.varacc<.3 || rst.sON.varacc<.3 || rst.mOFF.varacc<.3 || rst.mON.varacc<.3
        
       rst.sOFF.param(p,:) = rst.sOFF.param(p,:)*NaN; 
       rst.sON.param(p,:) = rst.sON.param(p,:)*NaN; 
       rst.mOFF.param(p,:) = rst.sOFF.param(p,:)*NaN; 
       rst.mON.param(p,:) = rst.sON.param(p,:)*NaN;        
       rst.Env.param(p,:) = rst.Env.param(p,:)*NaN; 
       rst.EnvON.param(p,:) = rst.EnvON.param(p,:)*NaN; 
       rst.EnvOFF.param(p,:) = rst.EnvOFF.param(p,:)*NaN; 
       rst.EnvM.param(p,:) = rst.EnvM.param(p,:)*NaN; 
       rst.EnvS.param(p,:) = rst.EnvS.param(p,:)*NaN; 
       
       rst.sOFF.ffit(p,:) = rst.sOFF.ffit(p,:)*NaN; 
       rst.sON.ffit(p,:) = rst.sON.ffit(p,:)*NaN; 
       rst.mOFF.ffit(p,:) = rst.sOFF.ffit(p,:)*NaN; 
       rst.mON.ffit(p,:) = rst.sON.ffit(p,:)*NaN; 
       rst.Env.ffit(p,:) = rst.Env.ffit(p,:)*NaN; 
       rst.EnvON.ffit(p,:) = rst.EnvON.ffit(p,:)*NaN; 
       rst.EnvOFF.ffit(p,:) = rst.EnvOFF.ffit(p,:)*NaN; 
       rst.EnvM.ffit(p,:) = rst.EnvM.ffit(p,:)*NaN; 
       rst.EnvS.ffit(p,:) = rst.EnvS.ffit(p,:)*NaN; 
       
%        Emat(p,:) = Emat(p,:)*NaN; 
%        Emat2(p,:) = Emat2(p,:)*NaN; 

    end
end
rst.sOFF.param(:,1) = rst.sOFF.param(:,1) + posdom(1);
rst.sON.param(:,1) = rst.sON.param(:,1) + posdom(1);
rst.mOFF.param(:,1) = rst.mOFF.param(:,1) + posdom(1);
rst.mON.param(:,1) = rst.mON.param(:,1) + posdom(1);

rst.Env.param(:,1) = rst.Env.param(:,1) + posdom(1);
rst.EnvON.param(:,1) = rst.EnvON.param(:,1) + posdom(1);
rst.EnvOFF.param(:,1) = rst.EnvOFF.param(:,1) + posdom(1);
rst.EnvS.param(:,1) = rst.EnvS.param(:,1) + posdom(1);
rst.EnvS.param(:,1) = rst.EnvS.param(:,1) + posdom(1); 
% 
% figure, hist(paramS(:,1)-paramM(:,1))
% xlabel('Spos - Mpos (degrees)')


%%

[dum posid] = sort(rst.Env.param(:,1)); %sort by retinotopic position
Ngood = length(find(~isnan(rst.Env.param(:,1)))); %Only plot good fits.



Mmap = zeros(64,3); 
scale = 0:63;
scale = scale.^.5 + .5;
scale = scale/max(scale)*63;
Mmap(:,2) = scale/63;

Smap = [64 0 64]/64; 
scale = 0:63;
scale = scale.^.5 ;
scale = scale/max(scale)*63;
Smap = scale'*Smap/63;

SMmap = (Smap)+flipud(Mmap);
SMmap = SMmap/max(SMmap(:));

SMmap = jet;
SMmap(:,1) = .4;
SMmap = SMmap(1:40,:);
SMmap = interp1((1:40)',SMmap,linspace(1,40,64)');
SMmapdum = SMmap;
SMmap(:,2) = SMmapdum(:,3);
SMmap(:,3) = SMmapdum(:,2);
%SMmap = flipud(SMmap);

clim = [0 1];

figure,
for pdum = 1:Ngood 
    p = posid(pdum);
    
    ax1 = subplot(1,Ngood+1,pdum);
   
    mdum = rst.EnvM.ffit(p,:);
    sdum = rst.EnvS.ffit(p,:);
    %mdum = mdum/max(mdum);
    %sdum = sdum/max(sdum);
    
    pS = sdum./(sdum+mdum);
    pS(find(pS>2 | pS<0)) = NaN;
    %plot(posdom,pS,'b'),ylim([0 1])
    %Win = mdum+sdum;
    Win = rst.Env.ffit(p,:)./rst.Env.param(p,3);
%     Win = Win-min(Win);
 %    Win = Win/max(Win);
    Win(find(Win<.37)) = 0;  %about 1.5sigma
    Win = sign(Win);
    %imagesc(1,posdom,pS','AlphaData',Win',[0 1]), colormap jet

    imagesc(1,posdom,pS','AlphaData',Win',[0 1]), axis image, %colorbar
    colormap(ax1,SMmap)
    caxis(clim)
    %caxis([0 1])
    axis xy
    
    
    ylim([-35 10])
    axis off
        
end

Win = nanmean(rst.Env.kern);
% Win = Win-min(Win);
Win = Win/max(Win);
Win(find(Win<.15)) = 0;
Win = sign(Win);
%Win = sign(Win);

%Make the composite: i.e. reconstruction
sN = rst.EnvS.ffit-rst.EnvS.param(:,4)*ones(1,length(posdom)); %subtract baseline from fit
mN = rst.EnvM.ffit-rst.EnvM.param(:,4)*ones(1,length(posdom));
sN = sN./(rst.Env.param(:,3)*ones(1,length(posdom))); %normalize them both by the amplitude of the envelope
mN = mN./(rst.Env.param(:,3)*ones(1,length(posdom)));

%Linear estimate at each retinotopic position:
good = find(~isnan(sum(rst.Env.param,2)));
clear pS_hat
for v = 1:length(posdom)
   H = sN(good,v)+mN(good,v);
   y = sN(good,v);
   
   H = sN(good,v);
   y = mN(good,v);

   pS_hat(v) = inv(H'*H)*H'*y;  
end

%Or just take the average of the ratio
% RetEnvelope = nanmean(sN + mN)/2;
% SEnvelope = nanmean(sN);
% MEnvelope = nanmean(mN);
% pS_hat = SEnvelope./(SEnvelope + MEnvelope);

ax1 = subplot(1,Ngood+1,pdum+1); imagesc(1,posdom,pS_hat','AlphaData',Win',[0 1]),
colormap(ax1,SMmap)
caxis(clim)

set(gca,'YaxisLocation','right')
ylim([-35 10])
axis xy
set(gca,'Xtick',[])

colorbar

figure,
for pdum = 1:Ngood
    p = posid(pdum);
    subplot(1,Ngood,pdum);

%     plot(rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:),posdom,'g'), hold on, 
%     hold on
%     plot(rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:),posdom,'b'),
%     dum = [rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:) rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:)];

    plot(rst.EnvON.ffit(p,:)-rst.EnvOFF.ffit(p,:),posdom,'k'), hold on, 
    dum = [rst.EnvON.ffit(p,:)-rst.EnvOFF.ffit(p,:)];


    xlim([min([dum(:); 0])-.001 max(dum(:))])
    hold on,
    plot([0 0], [posdom(1) posdom(end)])
    ylim([-30 5])
    axis off    
end
%%

figure,
for pdum = 1:Ngood
    p = posid(pdum);
    subplot(1,Ngood,pdum);

    plot(rst.mON.ffit(p,:),posdom,'g'), hold on, 
    hold on
    plot(-rst.sOFF.ffit(p,:),posdom,'b'),

%     dum = [rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:) rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:)];
%     xlim([min([dum(:); 0])-.001 max(dum(:))])
    axis off    
end
%%
id = find(Win);
[SWang] = getRetinaGradient(posdom) %Wang et al
figure, 
subplot(1,3,1), plot(posdom(id),SWang(id),'k'), ylim([10 100]), ylabel('%S Retinal ganglion cells')
subplot(1,3,3), plot(posdom(id),pS_hat(id)*100,'k'),  ylim([10 100]), ylabel('%S Cortical reconstruction')
subplot(1,3,2), semilogy(posdom(id),pS_hat(id)./SWang(id)*100,'k'), ylabel('S/M pooling bias')
xlabel('vertical retinotopy (deg)')

%%

vmat = ones(length(rst.Env.ffit(:,1)),1)*posdom;
[SWang] = getRetinaGradient(vmat)/100; %Wang et al

wS = rst.sON.ffit - rst.sON.param(:,4)*ones(1,size(SWang,2));
wS = wS./ (rst.sON.param(:,3)*ones(1,size(SWang,2)));
wS = wS./SWang;
wM = rst.mON.ffit - rst.mON.param(:,4)*ones(1,size(SWang,2));
wM = wM./ (rst.mON.param(:,3)*ones(1,size(SWang,2)));
wM = wM./(1-SWang);
% figure,
% subplot(1,2,1), plot(wS,wM,'.k')
% xlabel('SON weight'),ylabel('MON weight')

figure,
for pdum = 1:Ngood
    p = posid(pdum);
    subplot(10,8,pdum);

    plot(posdom,wS(p,:)/max(wS(p,:)),'b'), hold on, 
    hold on
    plot(posdom,wM(p,:)/max(wM(p,:)),'g')

%     dum = [rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:) rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:)];
%     xlim([min([dum(:); 0])-.001 max(dum(:))])
    axis off    
end

wS = rst.sOFF.ffit - rst.sOFF.param(:,4)*ones(1,size(SWang,2));
wS = wS./ (rst.sOFF.param(:,3)*ones(1,size(SWang,2)));
wS = wS./SWang;
wM = rst.mOFF.ffit - rst.mOFF.param(:,4)*ones(1,size(SWang,2));
wM = wM./ (rst.mOFF.param(:,3)*ones(1,size(SWang,2)));
wM = wM./(1-SWang);
% subplot(1,2,2), plot(wS,wM,'.k')
% xlabel('SOFF weight'),ylabel('MOFF weight')

figure,
for pdum = 1:Ngood
    p = posid(pdum);
    subplot(10,8,pdum);

    plot(posdom,wS(p,:)/max(wS(p,:)),'b'), hold on, 
    hold on
    plot(posdom,wM(p,:)/max(wM(p,:)),'g'),

%     dum = [rst.mON.ffit(p,:)-rst.mOFF.ffit(p,:) rst.sON.ffit(p,:)-rst.sOFF.ffit(p,:)];
%     xlim([min([dum(:); 0])-.001 max(dum(:))])
    axis off    
end


%%
vmat = ones(length(rst.Env.ffit(:,1)),1)*posdom;
[SWang] = getRetinaGradient(vmat)/100; %Wang et al

RSON = rst.sON.ffit - rst.sON.param(:,4)*ones(1,size(SWang,2));
RSON = RSON./ (rst.sON.param(:,3)*ones(1,size(SWang,2)));
RMON = rst.mON.ffit - rst.mON.param(:,4)*ones(1,size(SWang,2));
RMON = RMON./ (rst.mON.param(:,3)*ones(1,size(SWang,2)));

figure
subplot(3,2,1)
errorbar(posdom,nanmean(RMON),nanstd(RMON)/sqrt(size(RMON,1)),'g'), hold on, 
errorbar(posdom,nanmean(RSON),nanstd(RSON)/sqrt(size(RSON,1)),'b'), ylabel('R(v) = w(v)O(v) ON')
title('ON')

wSON = RSON./SWang;
wMON = RMON./(1-SWang);

ma = max(wMON,[],2)*ones(1,length(posdom));
wMON = wMON./ma;
ma = max(wSON,[],2)*ones(1,length(posdom));
wSON = wSON./ma;

subplot(3,2,5)
errorbar(posdom,nanmean(wMON),nanstd(wMON)/sqrt(size(wMON,1)),'g'), hold on, 
errorbar(posdom,nanmean(wSON),nanstd(wSON)/sqrt(size(wSON,1)),'b'), ylabel('w(v) ON')

RSOFF = rst.sOFF.ffit - rst.sOFF.param(:,4)*ones(1,size(SWang,2));
RSOFF = RSOFF./ (rst.sOFF.param(:,3)*ones(1,size(SWang,2)));
RMOFF = rst.mOFF.ffit - rst.mOFF.param(:,4)*ones(1,size(SWang,2));
RMOFF = RMOFF./ (rst.mOFF.param(:,3)*ones(1,size(SWang,2)));

subplot(3,2,2)
errorbar(posdom,nanmean(RMOFF),nanstd(RMOFF)/sqrt(size(RMOFF,1)),'g'), hold on, 
errorbar(posdom,nanmean(RSOFF),nanstd(RSOFF)/sqrt(size(RSOFF,1)),'b'), ylabel('R(v) = w(v)O(v) OFF')
title('OFF')

wSOFF = RSOFF./SWang;
wMOFF = RMOFF./(1-SWang);

ma = max(wMOFF,[],2)*ones(1,length(posdom));
wMOFF = wMOFF./ma;
ma = max(wSOFF,[],2)*ones(1,length(posdom));
wSOFF = wSOFF./ma;
subplot(3,2,6),
errorbar(posdom,nanmean(wMOFF),nanstd(wMOFF)/sqrt(size(wMOFF,1)),'g'), hold on, 
errorbar(posdom,nanmean(wSOFF),nanstd(wSOFF)/sqrt(size(wSOFF,1)),'b'), ylabel('w(v) OFF')


subplot(3,2,3)
plot(posdom,SWang,'b'), hold on,
plot(posdom,1-SWang,'g')
ylabel('RGC output; O(v)'),xlabel('vertical retinotopy')
subplot(3,2,4)
plot(posdom,SWang,'b'), hold on,
plot(posdom,1-SWang,'g')
ylabel('RGC output; O(v)'),xlabel('vertical retinotopy')


%%

colorVec = rst.mON.param(:,3)-rst.mOFF.param(:,3) + 1i*(rst.sON.param(:,3)-rst.sOFF.param(:,3));
hdom = [-180:10:180];
h = hist(angle(colorVec)*180/pi,hdom);
figure,bar(hdom,h)

%%
good = [3 4 5 6 11 14 15 23 27];
figure
for pid = 1:length(good)
    p = good(pid);
    Sprof = RSON(p,:)-RSOFF(p,:);
    Mprof = RMON(p,:)-RMOFF(p,:);
    
    %Sprof = wSON(p,:)-wSOFF(p,:);
    %Mprof = wMON(p,:)-wMOFF(p,:);
    
    %Mprof = Mprof/norm(Mprof);
    %Sprof = Sprof/norm(Sprof);
    
    proj = corrcoef(Sprof,Mprof);
    
    subplot(3,6,pid)
    %plot(Mprof,'g'), hold on, plot(Sprof,'b')

    title(num2str(proj(1,2)))
plot(RSON(p,:),'b'), hold on,
plot(-RSOFF(p,:),'b'), hold on,
plot(RMON(p,:),'g'), hold on,
plot(-RMOFF(p,:),'g')

% plot(wSON(p,:),'b'), hold on,
% plot(-wSOFF(p,:),'b'), hold on,
% plot(wMON(p,:),'g'), hold on,
% plot(-wMOFF(p,:),'g')
    
end



%% Plot vertical retinotopy
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
vrange = [posdom(1) posdom(end)];

vertplot = (rst.Env.param(:,1)-vrange(1))/(vrange(2)-vrange(1)) ;

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
title('vertical')
colorbar(ax1,'Ticks',[0 .5 1],'TickLabels',[posdom(1) 0 posdom(end)])
ylabel('microns')
%
% Fit a plane

vert_map = rst.Env.param(:,1);
cellLoc = CoM;

id = find(isnan(vert_map));
cellLoc(id,:) = [];
vert_map(id) = [];

[xmesh ymesh] = meshgrid(xdom,ydom);
H = [cellLoc(:,2)*xmicperpix cellLoc(:,1)*xmicperpix ones(size(cellLoc,1),1)];

vslope = inv(H'*H)*H'*vert_map;
vhatIm =  xmesh*vslope(1) + ymesh*vslope(2) + vslope(3);

subplot(2,2,2),imagesc(xdom,ydom,vhatIm,[posdom(1) posdom(end)]), colorbar, colormap jet, axis image
colorbar('Ticks',[posdom(1) 0 posdom(end)],'TickLabels',[posdom(1) 0 posdom(end)])
vmag = sqrt(sum(vslope(1:2).^2));
title(['Mag = ' num2str(round(1000*vmag)) 'deg/mm'])
xlabel('microns')

% Plot %S map

ES = sum(rst.EnvS.ffit,2);
EM = sum(rst.EnvM.ffit,2);
ES = rst.EnvS.param(:,3);
EM = rst.EnvM.param(:,3);
pS = ES./(EM+ES);
pS(find(pS>1 | pS<0)) = NaN;
pSidx = ceil(pS*64);
cdom = SMmap;

for p = 1:length(pS)
    
    if ~isnan(pS(p))
        
        vc = cdom(pSidx(p),:);
        
        ax2 = subplot(2,2,3);
        plot(CoM(p,2)*xmicperpix,CoM(p,1)*xmicperpix,'.','Color',vc,'MarkerSize',30)
        hold on
        xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis image
    end
    
end

xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)]), axis ij
title('%S')
vR = [0 1];
colorbar(ax2,'Ticks',[0 .5 1],'TickLabels',[0 .5 1])
colormap(ax2,cdom)

subplot(2,2,4),
scatter(rst.Env.param(:,1),pS)

id = find(~isnan(pS.*rst.Env.param(:,1)));
[r p] = corrcoef(pS(id),rst.Env.param(id,1))

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2))])
xlabel('RF position'), ylabel('%S')