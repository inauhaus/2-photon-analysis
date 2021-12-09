pF0

%initialize the Gui

global G_handles Analyzer cellS maskS idExamp

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

load('C:\2pScanboxAnalysis\ExpLogArray9','ExpLogArray')  

%A note on baselines:  I am not using the lowest two baselines [1/64 and
%1/128] because when I ran a simulation of the digital error (bit depth is
%7), there was a ton of cross contamination.  There was a big jump between
%1/64 and 1/32.  i.e. 1/32 was ok.  Also, exvec has a NaN if the experiment
%for that baseline does not exist... often the case for highest baseline
%with double current.

%%
clear animAll expt_mask exvec Kal basedom IsomRate idExampAll

%good
% idx = 1;
% animAll{idx} = 'se3'; %1/22/18
% expt_mask{idx} = 'u001_002'; %ori experiment
% %exvec{idx} = [NaN 9 8 7 6 5 4 NaN]; %could use 10 instead of 4
% exvec{idx} = [8 7 6 5 4 NaN]; %could use 10 instead of 4
% Kal{idx} = 'u001_011'; %retinotopy
% basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
% IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
% idExampAll{idx} = [];

%very few cells; Crappy retinotopy
idx = 1;
animAll{idx} = 'rl9'; %1/22/18
expt_mask{idx} = 'u002_003'; %ori experiment
%exvec{idx} = [NaN 10 9 7 6 5 4 NaN]; %Could use 11 instead of 4
exvec{idx} = [9 7 6 5 4 NaN]; %Could use 11 instead of 4
Kal{idx} = 'u002_002'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%sucks.  Crappy responses.
idx = idx+1;
animAll{idx} = 'rl7'; %1/22/18
expt_mask{idx} = 'u001_002'; %ori experiment
%exvec{idx} = [NaN 8 7 6 5 4 3 NaN];
exvec{idx} = [7 6 5 4 3 NaN];
Kal{idx} = 'u001_001'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%Really good.  Nice retinotopy and baseline trend
idx = idx+1;
animAll{idx} = 'rl7'; 
expt_mask{idx} = 'u001_013';
%exvec{idx} = [NaN 19 18 17 16 15 14 NaN]; 
exvec{idx} = [18 17 16 15 14 NaN]; %no baseline 6
Kal{idx} = 'u001_012'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%very few cells
idx = idx+1;
animAll{idx} = 'rl5'; 
expt_mask{idx} = 'u001_003'; %ori experiment
%exvec{idx} = [NaN 11 10 9 8 7 4 12]; %Could use 5 instead of 4 
exvec{idx} = [10 9 8 7 4 12]; %Could use 5 instead of 4 
Kal{idx} = 'u001_002'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

idx = idx+1;
animAll{idx} = 'rl7';
expt_mask{idx} = 'u002_003'; %ori experiment
%exvec{idx} = [NaN 16 14 12 10 9 5 6];  %Lots of other ones too
exvec{idx} = [14 12 10 9 5 6];  %Lots of other ones too
Kal{idx} = 'u002_021'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%Really good baseline trend.  All upper fields, so weak %S map
idx = idx+1;
animAll{idx} = 'rm2'; 
expt_mask{idx} = 'u000_002'; %ori experiment
%exvec{idx} = [13 11 10 9 8 7 6 4]; %could use 5 instead of 4
exvec{idx} = [10 9 8 7 6 4]; %could use 5 instead of 4
Kal{idx} = 'u000_001'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%Really good baseline and trend and %s gradient.  
idx = idx+1;
animAll{idx} = 'rm3'; %1/22/18
expt_mask{idx} = 'u000_002'; %ori experiment
%exvec{idx} = [17 15 14 12 10 8 6 3];  %could use 7 instead of 6.
exvec{idx} = [14 12 10 8 6 3];  %could use 7 instead of 6.
Kal{idx} = 'u000_001'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
%idExampAll{idx} = [5 10 25 30 55 60 100];
idExampAll{idx} = [5 10 30 100];

%good.  But lowest baseline creeps up in %S and has with the reverse gradient.
idx = idx+1;
animAll{idx} = 'rm4'; %1/23/18
expt_mask{idx} = 'u001_006';  %ori experiment
%exvec{idx} = [NaN 15 14 13 12 11 9 7]; %others to shose
exvec{idx} = [14 13 12 11 9 7]; %others to shose
Kal{idx} = 'u001_005'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

% idx = idx+1;
% animAll{idx} = 'rm6';
% expt_mask{idx} = 'u001_003';  %ori experiment
% exvec{idx} = [12 11 10 9 8 7 6 5];
% exvec{idx} = [10 9 8 7 6 5];
% Kal{idx} = 'u001_001'; %retinotopy
% basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
% IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
% idExampAll{idx} = [];

%Crappy retinotopy. Ok baseline trend... goes up at the end
idx = idx+1;
animAll{idx} = 'rm3'; %1/22/18
expt_mask{idx} = 'u001_004'; %ori experiment
%exvec{idx} = [12 11 10 9 8 6 5 4];  %could use 7 instead of 6.
exvec{idx} = [10 9 8 6 5 4];
Kal{idx} = 'u001_002'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%good trends, but %S is higher at lower baseline
idx = idx+1;
animAll{idx} = 'rm4'; %1/23/18
expt_mask{idx} = 'u002_005';  %ori experiment
%exvec{idx} = [12 11 10 9 8 7 6 5];
exvec{idx} = [10 9 8 7 6 5];
Kal{idx} = 'u002_002'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6);; %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
idExampAll{idx} = [];

%Really nice retinotopy.  Issac ran wrong color directions for lowest two
%baselines in experiments 16 and 15, so I have removed those.
idx = idx+1;
animAll{idx} = 'rm5'; 
expt_mask{idx} = 'u002_002'; %ori experiment
exvec{idx} = [11 10 9 8 5 4]; 
Kal{idx} = 'u002_003'; %retinotopy
basedom{idx} = logspace(log10(1/32),log10(1),6); %Highest has current values  UV250; g244
IsomRate{idx} = round(getIsomerizationRate(animAll{idx},Kal{idx}(4),2))/1000;
%idExampAll{idx} = [5 20 25 30 95 105];
idExampAll{idx} = [5 25 30 105];


% 

%% Rodless

clear animAll expt_mask exvec Kal basedom IsomRate idExampAll

idx = 1;
%idx = idx+1
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

exdom = 1:11;
exdom = 7
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
    
    %Loop each baseline
    Nb = length(Cexpt{ex});
    for bL = 1:Nb
        
        expt = Cexpt{ex}{bL}; %color experiment
        if ~isempty(expt)
            traceroot = 'C:\CellTraces\';
            tracepath = [traceroot anim '_' expt '_cellS'];
            try
                load([tracepath '_aligned to ' expt_mask{ex}(2:end)])
            catch
                load(tracepath,'cellS')
            end
            
            
            set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
            %             try
            %                 dataRoot = 'e:\2p_data\';
            %                 set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
            %                 Gsetdirectories %Load raw experiment and analyzer file
            %             catch
            %                 try
            %                     dataRoot = 'h:\2p_data\';
            %                     set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
            %                     Gsetdirectories %Load raw experiment and analyzer file
            %                 catch
            %                     dataRoot = 'g:\2p_data\';
            %                     set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
            %                     Gsetdirectories %Load raw experiment and analyzer file
            %
            %                 end
            %             end
            
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
            
            rst{ex}{bL} = ColorOri(1,2);
            
            pM{ex}{bL} = rst{ex}{bL}.pM;  %easier to have it organized like this
            M_Resp{ex}{bL} = rst{ex}{bL}.M_Resp;
            S_Resp{ex}{bL} = rst{ex}{bL}.S_Resp;
        else
            rst{ex}{bL} = [];
            pM{ex}{bL} = [];
            M_Resp{ex}{bL} = [];
            S_Resp{ex}{bL} = [];
        end
        
        
    end
   
    plotBaselineVColorExamples(rst{ex},IsomRate{ex},idExampAll{ex})
   
    
    
    %%%%%%%%%%%%%%%%%Now get Kalatsky info and compare to S/M %%%%%%%%%%%%%%%%
    expt = Kal{ex}; %Kalatsky
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS_aligned to ' expt_mask{ex}(2:end)];
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
    
    
    plot_pMvBaseRet_2photon(pM{ex},Kret{ex}.kmap_hor,Kret{ex}.kmap_vert,idExampAll{ex})
end


%%


%% Get the binning edges to use for both WT and gnat

edom = 1:11; %WT
edom = 1:3
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


prcdom(6) = 82;
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
edom = 1:3
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

id = find(pMAll<-.5 | pMAll>1.5);
pMAll(id) = NaN;
id = find(vretAll<-50 | vretAll>50);
vretAll(id) = NaN;

retBinEdge = [-50 0 30 60]; %Bin edges for retinotopy
clear muret mupM SEret SEpM mu_M mu_S SE_M SE_S

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

%% Data yield
edom = 1:3;
Ngood = zeros(1,6);
Nbad = zeros(1,6);
for eid = 1:length(edom) %loop each animal
    e = edom(eid);
    for bL = 1:length(pM{e}) %loop each baseline
        
        if ~isempty(pM{e}{bL})
            
            idBad = find(isnan(pM{e}{bL}.*Kret{e}.kmap_vert));
            idGood = find(~isnan(pM{e}{bL}.*Kret{e}.kmap_vert));
            
            Ngood(bL) = Ngood(bL)+length(idGood);
            Nbad(bL) = Nbad(bL)+length(idBad);
        end
        
    end
end

Ngood
Ntotal = Nbad+Ngood
Ngood./Ntotal
% 
% Ntotal = length(MRespAll);
% idBadGratings = find(isnan(MRespAll));
% idBadRetinotopy = find(isnan(vretAll));
% Bads = unique([idBadGratings idBadRetinotopy]);
% percYield = 1-length(Bads)/Ntotal



%% Plot %M vs. baseline, and fit the model


%subplot(1,length(prcdom)-1+1,bL+1),
figure(98)
for r = 1:length(vertdom)

    %plot(basedomPlot,mupM(:,r),'o','Color',ret2color(r,:),'MarkerSize',10),
    
    %errorbar(basedomPlot,(mupM(:,r)),SEpM(:,r),'.','Color',ret2color(r,:),'MarkerSize',10)
    errorbar(basedomPlot,(mupM(:,r)),SEpM(:,r),'.','Color',[0 0 0],'MarkerSize',10)
    hold on
    
    H = [basedomPlot' ones(length(basedomPlot),1)];
    y = log(mupM(:,r));
    xhat = inv(H'*H)*H'*y(:);  %xhat(1)*base + xhat(2);  exp(xhat(2))*exp(xhat(1)*base)
    
    [param ffit varaccount] = Expfit4(basedomPlot',mupM(:,r));
    %(1-param(1))*exp(-param(2)*domu) + param(1);
    
    
    baseModel = linspace(basedomPlot(1),basedomPlot(end),50);
    pMModel = exp(-param(2)*(baseModel).^param(1));
    
    
    
    %plot(basedum,exp(basedum*xhat(1)+xhat(2)),'Color',ret2color(r,:))
    plot(baseModel,pMModel,'Color',[0 0 0])
    hold on,
    plot([0 basedomPlot(end)],[1 1],'--k')
    hold on,
    plot([0 basedomPlot(end)],[0 0],'--k')
    
end
%set(gca,'ytick',[0 1])
xlabel('baseline (R*/rod/sec) x10^3')
ylabel('%M-opsin')
ylim([-.5 1.5])
set(gca,'XTick',round(basedomPlot))
xlim([0 400])


%%
vertIDall = [];
pSIDall = [];

%colorMat = [.025 .05 .1 .2 .35 .7]'*ones(1,3);
colorMat = [0 0 0 0 0 0]'*ones(1,3);
figure(99)
subplot(1,2,1)
for bL = 1:length(vertID)
   
    pSdum = 1-pMID{bL}(:); %rawdata
    %pSdum = pSdum/(1-mupM(bL,3));
    
    [dum bid] = min(abs(basedomPlot(bL)-baseModel));
    pSdum = pSdum/(1-pMModel(bid));  %use the model fit
    
    vertIDall = [vertIDall; vertID{bL}(:)];
    pSIDall = [pSIDall; pSdum];
    
    plot(vertID{bL}(:),1-pSdum,'.','MarkerEdgeColor',colorMat(bL,:))
    hold on
    
end

[dum ids] = sort(vertIDall);
[param ffit varaccount] = Sigfit2(vertIDall(ids),1-pSIDall(ids),[10 1]);

MRetModelDomain = -50:50;
d = MRetModelDomain-param(1);
MRetModel = 1./(1 + exp(-d*param(2)));

plot(MRetModelDomain,MRetModel,'r','LineWidth',2)

ylim([-.5 1.5]),xlim([-50 50])
hold on,
plot([-45 45],[1 1],'--k')
hold on,
plot([-45 45],[0 0],'--k')

set(gca,'ytick',[0 1])

set(gca,'XTick',[-40 -20 0 20 40])
%%
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
    
    %plot(vertIDall(ids),1-ffit*(1-mupM(bL,3)),'r','LineWidth',2)
    %[dum bid] = min(abs(basedomPlot(bL)-baseModel));
    %pSdum = pSdum/(1-pMModel(bid));  %use the model fit
    
    [dum bid] = min(abs(basedomPlot(bL)-baseModel));
    percRodModelSample = pMModel(bid);
    %pSdum = pSdum/(1-pMModel(bid));  %use the model fit
    
    plot(MRetModelDomain,1-(1-MRetModel)*(1-pMModel(bid)),'r','LineWidth',2)
    
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


%% Compare model to Gnat data, and Wang Demb model

load('C:\Users\Ian Nauhaus\Desktop\Mouse rod saturation paper\gnat_dat.mat')
%gnat_head
%percS_gnat_re
%v_gnat

vretAll_G = [];
pMAll_G = [];

for i = 1:length(percS_gnat_re)-1
   pMAll_G = [pMAll_G percS_gnat_re{i}(:,5)']; 
   vretAll_G = [vretAll_G v_gnat{i} + gnat_head(i)];
end
clear muret_G mupM_G SEret_G SEpM_G
figure
for v = 1:length(retBinEdge)-1
    idretbin = find(vretAll_G>retBinEdge(v) & vretAll_G<retBinEdge(v+1));
    retdum = vretAll_G((idretbin));
    pMAlldum = pMAll_G((idretbin));
    muret_G(v) = median(retdum);
    mupM_G(v) = median((pMAlldum));
    SEret_G(v) = std(retdum);%/sqrt(length(retdum));
    SEpM_G(v) = std((pMAlldum))/sqrt(length(pMAlldum));
    
    plot(muret_G(v),mupM_G(v),'.r','MarkerSize',20)
    hold on
    plot([muret_G(v) muret_G(v)], [mupM_G(v)-SEpM_G(v)/2 mupM_G(v)+SEpM_G(v)/2],'r','LineWidth',2)
    
end

figure(98)
for i = 1:length(mupM_G)
    hold on
    plot([10 400],[mupM_G(i) mupM_G(i)])
end



figure(99),
subplot(1,2,2)
scatter(vretAll_G,percSAll_G,'.k')

ylim([-.5 1.5]),xlim([-50 50])
hold on,
plot([-45 45],[1 1],'--k')
hold on,
plot([-45 45],[0 0],'--k')

hold on
plot(MRetModelDomain,MRetModel,'r','LineWidth',2)

opticDiskOffset = 30;
vdom = -60:50
[Sper] = getRetinaGradient(vdom-opticDiskOffset);

hold on
plot(vdom,(1-Sper/100),'b')

ylabel('%M')
xlabel('vertical retinotopy')

    set(gca,'ytick',[0 1])
    
    set(gca,'XTick',[-40 -20 0 20 40])



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