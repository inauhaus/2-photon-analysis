pF0

%initialize the Gui

global G_handles Analyzer cellS maskS idExamp

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'f:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

load('C:\2pScanboxAnalysis\ExpLogArray10','ExpLogArray')  

%A note on baselines:  I am not using the lowest two baselines [1/64 and
%1/128] because when I ran a simulation of the digital error (bit depth is
%7), there was a ton of cross contamination.  There was a big jump between
%1/64 and 1/32.  i.e. 1/32 was ok.  Also, exvec has a NaN if the experiment
%for that baseline does not exist... often the case for highest baseline
%with double current.

%%

%% Rodless

%clear animAll expt_mask exvec Kal basedom IsomRate idExampAll

%idx = 1;
idx = idx+1
%rodless
animAll{idx} = 'rn1'; %this seems to be the best of the rodless
expt_mask{idx} = 'u001_003'; %ori experiment
exvec{idx} = [12 11 10 9 7 5];  %Other options for last two baselines
Kal{idx} = 'u001_020'; %retinotopy
%basedom{idx} = logspace(log10(.0078125),log10(1),8); %Highest has current values  UV250; g244
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'rn2';   %Rodless mouse
expt_mask{idx} = 'u000_002';
exvec{idx} = [17 14 13 12 10 8 6 3]; %Other option for every baseline
exvec{idx} = [13 12 10 8 6 3]; %Other option for every baseline
Kal{idx} = 'u000_001'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [17 22];

idx = idx+1;
animAll{idx} = 'rn4';   %Rodless mouse
expt_mask{idx} = 'u001_007';
%exvec{idx} = [17 14 13 12 10 8 6 3]; %Other option for every baseline
exvec{idx} = [12 11 10 9 8 7]; %Other option for every baseline
Kal{idx} = 'u001_005'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
idExampAll{idx} = [7 12 22];
IsomRate{idx} = round(getIsomerizationRate('rm5',2,2))/1000;
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;


idx = 1;
animAll{idx} = 'rn1'; %5/16/19
expt_mask{idx} = 'u002_006'; %ori experiment  There is also an ori experiment at '1'
Cexpt{idx} = 'u002_006';  %S and M direction
Kal{idx} = 'u002_002'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'rn2'; %5/16/19
expt_mask{idx} = 'u001_008'; %ori experiment
Cexpt{idx} = 'u001_014';  %S and M direction; last experiment of day
Kal{idx} = 'u001_009'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'rn4'; %5/16/19
expt_mask{idx} = 'u003_003'; %S/M experiment
Cexpt{idx} = 'u003_003';
Kal{idx} = 'u003_002'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];


%%
clear animAll expt_mask exvec Kal basedom IsomRate idExampAll Cexpt

idx = 1;
animAll{idx} = 'se0'; %5/16/19
expt_mask{idx} = 'u002_006'; %ori experiment  There is also an ori experiment at '1'
Cexpt{idx} = 'u002_006';  %S and M direction
Kal{idx} = 'u002_002'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'ti3'; %5/16/19
expt_mask{idx} = 'u001_008'; %ori experiment
Cexpt{idx} = 'u001_014';  %S and M direction; last experiment of day
Kal{idx} = 'u001_009'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'ti4'; %5/16/19
expt_mask{idx} = 'u003_003'; %S/M experiment
Cexpt{idx} = 'u003_003';
Kal{idx} = 'u003_002'; %retinotopy
%IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];


%%
clear Cexpt
for i = 1:length(animAll)
    i
    for j = 1:length(exvec{i})
        if ~isnan(exvec{i}(j))
            Cexpt{i}{j} = ['u00' num2str(expt_mask{i}(4)) '_' sprintf('%03d',exvec{i}(j))];
        else
            Cexpt{i}{j} = [];
        end
    end
end


%% Load color x sf x ori expt
%clear pM rst Kret

global cellS maskS idExamp ACQinfo

exdom = 1:1;
%exdom = [7 11]; %use for paper figure 1
%exdom = [4 7 8 12]
%exdom = 1:3
clear rst pM
for eid = 1:length(exdom)
    
    ex = exdom(eid);
    
    idExamp = idExampAll{ex};
    
    anim = animAll{ex};
    
    expt = expt_mask{ex};
    try
        maskroot = 'C:\CellMasks\';
        maskpath = [maskroot anim '_' expt(1:8)];
        load(maskpath,'maskS')
    catch
        maskroot = 'C:\CellMasks_Issac\';
        maskpath = [maskroot anim '_' expt(1:8)];
        load(maskpath,'maskS')
    end
    
    expt = Cexpt{ex}; %color experiment
    
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    try
        load([tracepath '_aligned to ' expt_mask{ex}(2:end)])
    catch
        load(tracepath,'cellS')
    end
       
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    
    %Set dummy values for ACQinfo so that we don't need to have the
    %raw data files on hand
    ACQinfo.linesPerFrame = size(maskS.im{1},1);
    ACQinfo.pixelsPerLine = 676;
    ACQinfo.unblanked = 54:729;
    ACQinfo.msPerLine = 0.1261;
    ACQinfo.SBInfo.config.magnification = 1;
    %%%%%%%%%%%
    
    
    loadAnalyzer
    getCellStats
    
    rst{ex} = ColorOri(1,2);
    
    pM{ex} = rst{ex}.pM;  %easier to have it organized like this
    M_Resp{ex} = rst{ex}.M_Resp;
    S_Resp{ex} = rst{ex}.S_Resp;
    
end

  
    %%%%%%%%%%%%%%%%%Now get Kalatsky info and compare to S/M %%%%%%%%%%%%%%%%
    expt = Kal{ex}; %Kalatsky
    traceroot = 'C:\CellTraces\';
    %tracepath = [traceroot anim '_' expt '_cellS_aligned to ' expt_mask{ex}(2:end)];
    tracepath = [traceroot anim '_' expt '_cellS'];
    load(tracepath,'cellS')
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    %Gsetdirectories
    loadAnalyzer
    
    getCellStats
    
    %f1 = f1meanimage;  %Build F1 images (takes the longest)
    %First set dirs and hit 'Compute F0 images' within pF0
    f1 = CondF1_cellMask2(.9);    
    [kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing
    
    tilt = getTilt(anim,ExpLogArray);
    kmap_vert = kmap_vert+tilt;  
    
    projectorAdjustment = 1;
    
    if projectorAdjustment
        kmap_horx = kmap_hor;
        kmap_vertx = kmap_vert;
        kmap_hor = kmap_vertx;
        kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
    end
    
    plotMaskedRetinotopy(kmap_hor(:),kmap_vert(:))
    
    Kret{ex}.kmap_hor = kmap_hor;
    Kret{ex}.kmap_vert = kmap_vert;
    
    pMdum{1} = pM{ex};
    plot_pMvBaseRet_2photon(pMdum,Kret{ex}.kmap_hor,Kret{ex}.kmap_vert,idExampAll{ex})
end


%%


%% Get the binning edges to use for both WT and gnat

edom = 1:11; %WT
%edom = 12:14; %gnat
baseAll = [];
for eid = 1:length(edom)
    e = edom(eid);
    for bL = 1:length(pM{e})
        if ~isempty(pM{e}{bL})            
            
            %Calculation below is important: Index of 3 is for rods (as opposed to S and M-opsin)
            %Calibration was run at baseline 5
            baseAll = [baseAll basedom{1}(bL)*ones(1,length(pM{e}{bL}))*IsomRate{e}(3)*2]; %x2 because isom rate was basedom = 0.5
            %baseAll = [baseAll basedom{1}(bL)*ones(1,length(pM{e}{bL}))];
        end
    end
end

prcdom = linspace(0,100,7); %Bin edges for the baselines
%Loop each bin of light level and get stats.
for bL = 1:length(prcdom)
    
    baseWinEdgedum = prctile(log2(baseAll),prcdom(bL));
    %I was previously leaving out data points on an edge w/o doing this:
    if bL == 1
        baseWinEdgedum = 0;
    elseif bL == length(prcdom)
        baseWinEdgedum = 100;
    end
    baseWinEdge(bL) = baseWinEdgedum;
end

%% get the mean
clear basedomPlot
for bL = 1:length(prcdom)-1
    
    baseWin = [baseWinEdge(bL) baseWinEdge(bL+1)];  %baseWinEdge computed in cell above

    id = find(log2(baseAll) > baseWin(1) & log2(baseAll) <= baseWin(2));
    %id2 = find(~isnan(vretAll(id).*pMAll(id)));
    %id = id(id2);
    
    basedomPlot(bL) = geomean(baseAll(id));
end
    
   
%%

edom = 1:11; %WT
%edom = 12:14; %gnat

%Accumulate across experiments
pMAll = []; vretAll = []; hretAll = []; baseAll = []; MRespAll = []; SRespAll = [];
for eid = 1:length(edom)
    e = edom(eid);
    for bL = 1:length(pM{e})
        if ~isempty(pM{e}{bL})            
            pMAll = [pMAll pM{e}{bL}];
            MRespAll = [MRespAll M_Resp{e}{bL}];
            SRespAll = [SRespAll S_Resp{e}{bL}];
            vretAll = [vretAll Kret{e}.kmap_vert];
            hretAll = [hretAll Kret{e}.kmap_hor];
            
            %Calculation below is important: Index of 3 is for rods (as opposed to S and M-opsin)
            %Calibration was run at baseline 5
            baseAll = [baseAll basedom{1}(bL)*ones(1,length(pM{e}{bL}))*IsomRate{e}(3)*2]; %x2 because isom rate was basedom = 0.5
            %baseAll = [baseAll basedom{1}(bL)*ones(1,length(pM{e}{bL}))];
        end
    end
end

figure

prcdom = linspace(0,100,7); %Bin edges for the baselines

id = find(pMAll<-.5 | pMAll>1.5);
pMAll(id) = NaN;
id = find(vretAll<-50 | vretAll>50);
vretAll(id) = NaN;

retBinEdge = [-30 0 30 60]; %Bin edges for retinotopy
clear muret mupM SEret SEpM

%Loop each bin of light level and get stats.
for bL = 1:length(prcdom)-1
    
%     baseWin = [prctile(log2(baseAll),prcdom(bL)) prctile(log2(baseAll),prcdom(bL+1))]
%     %I was previously leaving out data points on an edge w/o doing this:
%     if bL == 1
%         baseWin(1) = 0;
%     elseif bL == 6
%         baseWin(2) = 100;
%     end

    baseWin = [baseWinEdge(bL) baseWinEdge(bL+1)];  %baseWinEdge computed in cell above

    id = find(log2(baseAll) > baseWin(1) & log2(baseAll) <= baseWin(2));
    id2 = find(~isnan(vretAll(id).*pMAll(id)));
    id = id(id2);
    
     %basedomPlot(bL) = geomean(baseAll(id));  
    
    muM(bL) = median(MRespAll(id));
    muS(bL) = median(SRespAll(id));
    sigM(bL) = std(MRespAll(id))/sqrt(length(id));
    sigS(bL) = std(SRespAll(id))/sqrt(length(id));
    
    vertID{bL} = vretAll(id);
    pMID{bL} = pMAll(id);
    
    N = length(id);
    subplot(1,length(prcdom)-1+1,bL),
    scatter(vretAll(id),pMAll(id),'.k'), 

    ylim([-.5 1.5]),xlim([-50 50])       
    hold on,
    plot([-45 45],[1 1],'--k')
    hold on,
    plot([-45 45],[0 0],'--k')
    
    set(gca,'ytick',[0 1])
    
    if bL ~= length(prcdom)-1
        set(gca,'xticklabels',[])
    end
    
    if bL == length(prcdom)-1
        [dum ids] = sort(vretAll(id));
        [param ffit varaccount] = Sigfit2(vretAll(id(ids)),pMAll(id(ids)),[15 3]);
        hold on
        vretfit = vretAll(id(ids));
        plot(vretfit,ffit,'r','LineWidth',2)
    end
    
    for v = 1:length(retBinEdge)-1
       idretbin = find(vretAll(id)>retBinEdge(v) & vretAll(id)<retBinEdge(v+1));
       retdum = vretAll(id(idretbin));
       pMAlldum = pMAll(id(idretbin));
       muret(bL,v) = median(retdum);
       mupM(bL,v) = median((pMAlldum));
       SEret(bL,v) = std(retdum);%/sqrt(length(retdum));
       SEpM(bL,v) = std((pMAlldum))/sqrt(length(pMAlldum));
       
       MAlldum = MRespAll(id(idretbin));
       SAlldum = SRespAll(id(idretbin));
       mu_M(bL,v) = mean((MAlldum));
       mu_S(bL,v) = mean((SAlldum));
       SE_M(bL,v) = std((MAlldum))/sqrt(length(MAlldum));
       SE_S(bL,v) = std((SAlldum))/sqrt(length(SAlldum));
             
       plot(muret(bL,v),mupM(bL,v),'.r','MarkerSize',20)
       hold on
       plot([muret(bL,v) muret(bL,v)], [mupM(bL,v)-SEpM(bL,v)/2 mupM(bL,v)+SEpM(bL,v)/2],'r','LineWidth',2)
       
    end
    
    hold on,

    %plot(muret(bL,:),mupM(bL,:),'.r','MarkerSize',20)
    
    isomRate = 2.^mean(log2(baseAll(id)));
    
    %title(['R*/rod/sec = ' num2str(round(isomRate)) ' x10^3']);
    xlim([-45 45])
    
    
end

retCdom = jet;
RetLims = [-50 50];

vertdom = (muret(end,:));

retCLUT = linspace(RetLims(1),RetLims(2),64);

%Make retinotoy color table
clear ret2color
%vertdomdum = vertdom;
vertdomdum = linspace(RetLims(1),RetLims(2),length(vertdom));
for q = 1:length(vertdom)
    [dum id] = min(abs(vertdomdum(q)-retCLUT));
    ret2color(q,:) = retCdom(id,:);
end
%% Plot %M vs. baseline.


%subplot(1,length(prcdom)-1+1,bL+1),
figure
for r = 1:length(vertdom)
    %plot(basedomPlot,mupM(:,r),'o','Color',ret2color(r,:),'MarkerSize',10),
    
    %errorbar(basedomPlot,(mupM(:,r)),SEpM(:,r),'.','Color',ret2color(r,:),'MarkerSize',10)
    errorbar(basedomPlot,(mupM(:,r)),SEpM(:,r),'.','Color',[0 0 0],'MarkerSize',10)
    hold on
    
    H = [basedomPlot' ones(length(basedomPlot),1)];
    y = log2(mupM(:,r));
    xhat = inv(H'*H)*H'*y(:);
    
    basedum = linspace(basedomPlot(1),basedomPlot(end),50);
    
    %plot(basedum,2.^(basedum*xhat(1)+xhat(2)),'Color',ret2color(r,:))
    plot(basedum,2.^(basedum*xhat(1)+xhat(2)),'Color',[0 0 0])
    hold on,
    plot([0 basedomPlot(end)],[1 1],'--k')
    hold on,
    plot([0 basedomPlot(end)],[0 0],'--k')
    
end
set(gca,'ytick',[0 1])
xlabel('baseline (R*/rod/sec) x10^3')
ylabel('%M-opsin')
ylim([-.5 1.5])
set(gca,'XTick',round(basedomPlot))


%% Plot %M vs. baseline.  (Linear)


%subplot(1,length(prcdom)-1+1,bL+1),
figure
for r = 1:length(vertdom)
    %plot(basedomPlot,mupM(:,r),'o','Color',ret2color(r,:),'MarkerSize',10),
    
    errorbar(log2(basedomPlot),(mupM(:,r)),SEpM(:,r),'.','Color',ret2color(r,:),'MarkerSize',10)
    hold on
    
    H = [log2(basedomPlot') ones(length(basedomPlot),1)];
    y = (mupM(:,r));
    xhat = inv(H'*H)*H'*y(:);
    
    basedum = log2(linspace(basedomPlot(1),basedomPlot(end),50));
    
    
    plot(basedum,(basedum*xhat(1)+xhat(2)),'Color',ret2color(r,:))
    hold on,
    %plot([0 log2(basedomPlot(end))],[1 1],'--k')
    hold on,
    %plot([0 log2(basedomPlot(end))],[0 0],'--k')
    
   
end
set(gca,'ytick',[0 1])
xlabel('baseline (R*/rod/sec) x10^3')
ylabel('%M-opsin')
ylim([-.5 1.5])
set(gca,'XTick',log2(round(basedomPlot)))
%%
vertIDall = [];
pSIDall = [];

%colorMat = [.025 .05 .1 .2 .35 .7]'*ones(1,3);
colorMat = [0 0 0 0 0 0]'*ones(1,3);
figure
for bL = 1:length(vertID)
   
    pSdum = 1-pMID{bL}(:);
    pSdum = pSdum/(1-mupM(bL,3));
    
    vertIDall = [vertIDall; vertID{bL}(:)];
    pSIDall = [pSIDall; pSdum];
    
    plot(vertID{bL}(:),1-pSdum,'.','MarkerEdgeColor',colorMat(bL,:))
    hold on
    
end

[dum ids] = sort(vertIDall);
[param ffit varaccount] = Sigfit2(vertIDall(ids),pSIDall(ids),[10 1]);

plot(vertIDall(ids),1-ffit,'r','LineWidth',2)

ylim([-.5 1.5]),xlim([-50 50])
hold on,
plot([-45 45],[1 1],'--k')
hold on,
plot([-45 45],[0 0],'--k')

set(gca,'ytick',[0 1])

set(gca,'XTick',[-40 -20 0 20 40])

% bLfit = 6;
% wt = (1-mupM(:,3));
% wt = wt/wt(bLfit);
% [dum ids] = sort(vertID{bLfit});
% [param ffit varaccount] = Sigfit2(vertID{bLfit}(ids),pMID{bLfit}(ids),[10 .5])
% hold on
% vretfit = vertID{bLfit}(ids);
% plot(vretfit,ffit,'r','LineWidth',2)

figure
for bL = 1:length(vertID)
    
    subplot(1,6,bL),
    
    scatter(vertID{bL},pMID{bL},'.k'), 
    hold on
    
   % pSall = (1-pMID{bL})
    
    %ffitdum = (1-ffit)*wt(bL);
    %ffitdum = 1-ffitdum;
    
    plot(vertIDall(ids),1-ffit*(1-mupM(bL,3)),'r','LineWidth',2)
    
    plot(muret(bL,1),mupM(bL,1),'.b','MarkerSize',20)
    plot(muret(bL,2),mupM(bL,2),'.g','MarkerSize',20)
    plot(muret(bL,3),mupM(bL,3),'.r','MarkerSize',20)
    
    ylim([-.5 1.5]),xlim([-50 50])       
    hold on,
    plot([-45 45],[1 1],'--k')
    hold on,
    plot([-45 45],[0 0],'--k')
    
    set(gca,'ytick',[0 1])
    
    set(gca,'XTick',[-40 -20 0 20 40])
end

%%  Absolute df/f vs baseline

% figure,
% subplot(2,1,1)
% plot(log(baseAll),(MRespAll+SRespAll)/2,'.k')
% hold on, errorbar(log(basedomPlot),(muM+muS)/2,(sigM+sigS)/2,'r')
% ylabel('dF/F S-iso + M-iso')
% 
% subplot(2,1,2),
% %plot(log(baseAll),SRespAll,'.k')
% hold on, errorbar(log(basedomPlot),muS,sigS,'b')
% ylabel('dF/F S-isolating')
% 
% %plot(log(baseAll),MRespAll,'.k')
% hold on, errorbar(log(basedomPlot),muM,sigM,'g')
% ylabel('dF/F M-isolating')
% ylim([-.05 .3])

for i = 1:length(basedomPlot)
   
    xlabs{i} = num2str(round(basedomPlot(i)*10)/10);
    
end

figure,
subplot(3,1,1)
errorbar(log(basedomPlot),mu_M(:,1),SE_M(:,1),'g')
hold on, errorbar(log(basedomPlot),mu_S(:,1),SE_S(:,1),'b')
%hold on, plot([-.25 4], [0 0],'--k')
set(gca,'Xtick',log(basedomPlot),'XtickLabel',xlabs)
ylabel('dF/F lower fields')
ylim([-.05 .3])

subplot(3,1,2)
errorbar(log(basedomPlot),mu_M(:,2),SE_M(:,2),'g')
hold on, errorbar(log(basedomPlot),mu_S(:,2),SE_S(:,2),'b')
set(gca,'Xtick',log(basedomPlot),'XtickLabel',xlabs)
ylabel('dF/F middle fields')
ylim([-.05 .3])
%hold on, plot([-.25 4], [0 0],'--k')

subplot(3,1,3)
errorbar(log(basedomPlot),mu_M(:,3),SE_M(:,3),'g')
hold on, errorbar(log(basedomPlot),mu_S(:,3),SE_S(:,3),'b')
set(gca,'Xtick',log(basedomPlot),'XtickLabel',xlabs)
ylabel('dF/F upper fields')
ylim([-.05 .3])
xlabel('R*/sec  (x10^3)')
%hold on, plot([-.25 4], [0 0],'--k')

figure
plot(log(basedomPlot),mu_M(:,3)./(mu_S(:,3)+mu_M(:,3)),'g')
ylim([0 .8])
%% Wang Demb comparison
figure

clear muret mupM SEret SEpM

bL = length(prcdom)-1;
baseWin = [prctile(log2(baseAll),prcdom(bL)) prctile(log2(baseAll),prcdom(bL+1))];
id = find(log2(baseAll) >= baseWin(1) & log2(baseAll) < baseWin(2));
id2 = find(~isnan(vretAll(id).*pMAll(id)));
id = id(id2);

basedomPlot(bL) = geomean(baseAll(id));


scatter(vretAll(id),pMAll(id),'.k'),
ylim([-.5 1.5]),xlim([-50 50])
hold on,
plot([-45 45],[1 1],'--k')
hold on,
plot([-45 45],[0 0],'--k')

opticDiskOffset = 30;
vdom = -40:40
[Sper] = getRetinaGradient(vdom-opticDiskOffset);

hold on
plot(vdom,(1-Sper/100),'b')

hold on
plot(vretfit,ffit,'r')

ylabel('%M')
xlabel('vertical retinotopy')

%%



figure,
subplot(2,3,1)
for bL = 1:length(basedom)
    plot(vertdom,pMMat(:,bL),'.-','LineWidth',baseLineWid(bL),'MarkerSize',baseLineWid(bL)*6,'Color',[0 0 0]), hold on
    xlabel('vertical retinotopy'), ylabel('%S')
    hold on
end


% for bL = 1:length(pM{e})
%     id = find(baseAll == basedom{1}(bL));
%     subplot(1,length(pM{e}),bL),
%     scatter(vretAll(id),pMAll(id),'.k'), ylim([-2 2]),xlim([-45 45])
% end

    

