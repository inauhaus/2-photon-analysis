pF0

%initialize the Gui

global G_handles Analyzer cellS maskS idExamp

set(G_handles.epistart,'String','500');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3250'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

clear rst0_90 rst45_135 Kret

%%
clear animAll expt_mask expt_0v90 expt_45v135 expt_Kal idExampAll

id = 1;

animAll{id} = 'rk5';
expt_mask{id} = 'u001_010';
expt_0v90{id} = 'u001_010';
expt_45v135{id} = 'u001_011';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [3 21 75 81 92];

%Very few responsive cells
id = id+1;
animAll{id} = 'rk7';
expt_mask{id} = 'u001_006';
expt_0v90{id} = 'u001_006';
expt_45v135{id} = 'u001_007';
expt_Kal{id} = 'u001_000';
idExampAll{id} = [3 21];

id = id+1;
animAll{id} = 'rl0';
expt_mask{id} = 'u001_005';
expt_0v90{id} = 'u001_005';
expt_45v135{id} = 'u001_006';
expt_Kal{id} = 'u001_001';
idExampAll{id} = [19 21 25 49 58 78 98];

id = id+1;
animAll{id} = 'rl1';
expt_mask{id} = 'u002_006';
expt_0v90{id} = 'u002_006';
expt_45v135{id} = 'u002_007';
expt_Kal{id} = 'u002_003';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rk3';
expt_mask{id} = 'u002_009';
expt_0v90{id} = 'u002_010';
expt_45v135{id} = 'u002_011';
expt_Kal{id} = 'u002_004';
idExampAll{id} = [];

%very few cells
id = id+1;
animAll{id} = 'rl2';
expt_mask{id} = 'u000_005';
expt_0v90{id} = 'u000_006';
expt_45v135{id} = 'u000_007';
expt_Kal{id} = 'u000_001';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rk5';
expt_mask{id} = 'u002_002';
expt_0v90{id} = 'u002_004';
expt_45v135{id} = 'u002_005';
expt_Kal{id} = 'u002_001';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rl0';
expt_mask{id} = 'u002_006';
expt_0v90{id} = 'u002_004';
expt_45v135{id} = 'u002_006';
expt_Kal{id} = 'u002_001';
idExampAll{id} = [];



%% New ones

%Can't get any signal in color vs. lum experiment
% animAll{1} = 'rk9'; %Nice retinotopy, but color yields no response.  Very low yield.
% expt_mask{1} = 'u001_002';
% expt_0v90{1} = 'u001_007';
% expt_45v135{1} = 'u001_008';
% expt_Kal{1} = 'u001_001';
% idExampAll{1} = [];

id = id+1;
animAll{id} = 'rm0';  %ok.
expt_mask{id} = 'u000_007';
expt_0v90{id} = 'u000_008';
expt_45v135{id} = 'u000_007';
expt_Kal{id} = 'u000_005';
idExampAll{id} = [];

%
% animAll{1} = 'rl9';
% expt_mask{1} = 'u001_004';
% expt_0v90{1} = 'u001_008';
% expt_45v135{1} = 'u001_007';
% expt_Kal{1} = 'u001_005';

id = id+1;
animAll{id} = 'rm2';  %Weak but consistent
expt_mask{id} = 'u001_008';
expt_0v90{id} = 'u001_007';
expt_45v135{id} = 'u001_008';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rm5'; %consistent
expt_mask{id} = 'u001_007';
expt_0v90{id} = 'u001_006';
expt_45v135{id} = 'u001_007';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rr2';  %Very consistent, but some negative responses
expt_mask{id} = 'u001_007';
expt_0v90{id} = 'u001_006';
expt_45v135{id} = 'u001_007';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [];

id = id+1;
animAll{id} = 'rr9';  %Weak
expt_mask{id} = 'u001_005';
expt_0v90{id} = 'u001_004';
expt_45v135{id} = 'u001_005';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [];

%Very few cells
% id = id+1;
% animAll{id} = 'rr8';  %Not consistent
% expt_mask{id} = 'u001_004';
% expt_0v90{id} = 'u001_003';
% expt_45v135{id} = 'u001_004';
% expt_Kal{id} = 'u001_001';
% idExampAll{id} = [];

id = id+1;
animAll{id} = 'rw5';  %Jackpot
expt_mask{id} = 'u001_009';
expt_0v90{id} = 'u001_007';
expt_45v135{id} = 'u001_009';
expt_Kal{id} = 'u001_003';
idExampAll{id} = [];

% id = id+1;
% animAll{id} = 'rw9';  %Very few cells.  Almost not worth using
% expt_mask{id} = 'u001_006';  %looks like a lot of brain movement
% expt_0v90{id} = 'u001_005';
% expt_45v135{id} = 'u001_006';
% expt_Kal{id} = 'u001_002';
% idExampAll{id} = [];

% id = id+1;
% animAll{id} = 'rx0';  %Very few cells
% expt_mask{id} = 'u001_005';
% expt_0v90{id} = 'u001_004';
% expt_45v135{id} = 'u001_005';
% expt_Kal{id} = 'u001_002';
% idExampAll{id} = [];

%almost no cells
% id = id+1;
% animAll{id} = 'rw8'; 
% expt_mask{id} = 'u001_010';
% expt_0v90{id} = 'u001_009';
% expt_45v135{id} = 'u001_010';
% expt_Kal{id} = 'u001_007';
% idExampAll{id} = [];

id = id+1;
animAll{id} = 'rw5';  %good. 
expt_mask{id} = 'u002_005';
expt_0v90{id} = 'u002_004';
expt_45v135{id} = 'u002_005';
expt_Kal{id} = 'u002_002';
idExampAll{id} = [];

%Something strange with this one.
% id = id+1;
% animAll{id} = 'rx1';  
% expt_mask{id} = 'u002_005';
% expt_0v90{id} = 'u002_004';
% expt_45v135{id} = 'u002_005';
% expt_Kal{id} = 'u002_002';
% idExampAll{id} = [];

%% This one has all four color directions
id = id+1; 
animAll{id} = 'rk4';  %good
expt_mask{id} = 'u001_003'; 
expt_0v90{id} = 'u001_005';
expt_45v135{id} = 'u001_005';
expt_Kal{id} = 'u001_002';
idExampAll{id} = [];

%% Load color x sf x ori expt

global cellS maskS idExamp

Kalflag = 1;

exdom = 1:length(expt_mask);
exdom = 1
exdom = 1:10
%exdom = 9;
for eid = 1:length(exdom)
    
    ex = exdom(eid);
    
    idExamp = idExampAll{ex};
    
    anim = animAll{ex};
    
    %%%%%%First get 0 and 90%%%%%%%%%%%%
    
    expt = expt_mask{ex};
    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    load(maskpath,'maskS')
    
    expt = expt_0v90{ex}; %S vs M
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    try
        load([tracepath '_aligned to ' expt_mask{ex}(2:end)])
    catch
        load(tracepath,'cellS')
    end
    
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    try        
        dataRoot  = 'h:\2p_data\';
        set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
        Gsetdirectories %Load raw experiment and analyzer file
        
    catch
        try
            dataRoot  = 'g:\2p_data\';
            set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
            Gsetdirectories %Load raw experiment and analyzer file
        catch
            try
                dataRoot  = 'e:\2p_data\';
                set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
                Gsetdirectories %Load raw experiment and analyzer file
            catch
                
                loadAnalyzer %Just load analyzer file
            end

        end
    end
    
    %%%%%%%%%%%
    
    getCellStats
    
    rst0_90{ex} = ColorSfreq4(1,2);
    
    %%%%%%Now get 45 and 135%%%%%%%%
    %Use mask from 001_010
    expt = expt_45v135{ex}; %color vs lum
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    try
        load([tracepath '_aligned to ' expt_mask{ex}(2:end)])
    catch
        load(tracepath,'cellS')
    end
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    try
        Gsetdirectories %Load raw experiment and analyzer file
    catch
        loadAnalyzer %Just load analyzer file
    end
    
    getCellStats
    
    rst45_135{ex} = ColorSfreq4(1,2);
    
    plotBothColorDirs(rst0_90{ex},rst45_135{ex})
    
    
    if Kalflag
        %%%%%%%%%%%%%%%%%Now get Kalatsky info and compare to S/M %%%%%%%%%%%%%%%%
        expt = expt_Kal{ex}; %Kalatsky
        traceroot = 'C:\CellTraces\';
        tracepath = [traceroot anim '_' expt '_cellS_aligned to ' expt_mask{ex}(2:end)];
        load(tracepath,'cellS')
        set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
        set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
        try
            Gsetdirectories %Load raw experiment and analyzer file
        catch
            loadAnalyzer %Just load analyzer file
        end
        
        getCellStats
        
        %f1 = f1meanimage;  %Build F1 images (takes the longest)
        %First set dirs and hit 'Compute F0 images' within pF0
        f1 = CondF1_cellMask2(.9);
        
        [kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing
        
        projectorAdjustment = 1;
        
        if projectorAdjustment
            kmap_horx = kmap_hor;
            kmap_vertx = kmap_vert;
            kmap_hor = kmap_vertx;
            kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
        end
        
        plotMaskedRetinotopy(kmap_hor(:),kmap_vert(:))
        
        figure, scatter(kmap_vert,rst0_90{ex}.pS,'k')
        ylabel('%S')
        xlabel('Vertical retinotopy')
        ylim([0 1])
        
        id = find(~isnan(kmap_vert(:).*rst0_90{ex}.pS(:)));
        [r p] = corrcoef(kmap_vert(id),rst0_90{ex} .pS(id));
        r = r(1,2); p = p(1,2);
        
        title(['p = ' num2str(p)])
        
        Kret{ex}.kmap_hor = kmap_hor;
        Kret{ex}.kmap_vert = kmap_vert;
    end
end

%%


%% Load color x sf x ori expt
global cellS maskS

anim = 'rk4'; %This experiment has all 4 color direction together
eid = 16;
eid = 11
%ok; sf bias is in correct direction.  Consistently low-pass in color
%direction
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_003';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_005'; %S,S+M,M,S-M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 001_003'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
try
    Gsetdirectories %Load raw experiment and analyzer file
catch
    loadAnalyzer %Just load analyzer file
end

getCellStats

c1 = 1; c2 = 3;
rst0_90{eid} = ColorSfreq2(c1,c2);

c1 = 2; c2 = 4;
rst45_135{eid} = ColorSfreq2(c1,c2);

%%%%%%%%%%%%%%%%%Now get Kalatsky info and compare to S/M %%%%%%%%%%%%%%%%
expt = 'u001_002'; %Kalatsky
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 001_003'];
load(tracepath,'cellS')

dataRoot = 'E:\2p_data\'
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])

Gsetdirectories

getCellStats

%f1 = f1meanimage;  %Build F1 images (takes the longest)
%First set dirs and hit 'Compute F0 images' within pF0
f1 = CondF1_cellMask2(.9);

[kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing

projectorAdjustment = 1;

if projectorAdjustment    
    kmap_horx = kmap_hor;
    kmap_vertx = kmap_vert;
    kmap_hor = kmap_vertx;
    kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
end

plotMaskedRetinotopy(kmap_hor(:),kmap_vert(:)) 

figure, scatter(kmap_vert,rst0_90{eid}.pS,'k')
ylabel('%S')
xlabel('Vertical retinotopy')
ylim([0 1])

id = find(~isnan(kmap_vert(:).*rst0_90{eid}.pS(:)));
[r p] = corrcoef(kmap_vert(id),rst0_90{eid}.pS(id));
r = r(1,2); p = p(1,2);

title(['p = ' num2str(p)])

Kret{eid}.kmap_hor = kmap_hor;
Kret{eid}.kmap_vert = kmap_vert;


%% Population measures

%% Accumulate data from all experiments

vars = {'omag', 'opref' ,'dmag' ,'dpref', 'sfpref', 'sfhco', 'sfBW', 'sfBP','sfCoM'} 
kmap_vert = [];
kmap_hor = [];
pSAll = [];
pColorAll = [];

for i = 1:length(rst45_135)-1
    
    
    
    for j = 1:length(vars)
        
        if i == 1
            eval([vars{j} 'Color = [];'])
            eval([vars{j} 'Lum = [];'])
            eval([vars{j} 'SOpsin = [];'])
            eval([vars{j} 'MOpsin = [];'])
        end
        
        
        eval([vars{j} 'Color = [' vars{j} 'Color rst45_135{i}.'  vars{j} '{2}];'])
        eval([vars{j} 'Lum = [' vars{j} 'Lum rst45_135{i}.'  vars{j} '{1}];'])
        
        eval([vars{j} 'SOpsin = [' vars{j} 'SOpsin rst0_90{i}.'  vars{j} '{2}];'])
        eval([vars{j} 'MOpsin = [' vars{j} 'MOpsin rst0_90{i}.'  vars{j} '{1}];'])
        
    end
    %                 catch
    %
    %             eval([vars{j} 'All = [' vars{j} 'All rst45_135{i}.'  vars{j} '];'])
    %
    %         end
    
    
    
    kmap_hor = [kmap_hor Kret{i}.kmap_hor];
    kmap_vert = [kmap_vert Kret{i}.kmap_vert];
    
end 
    


%Get all the tuning curves and fits
sftcAll_45 = [];
sftcAll_135 = [];
sftcAll_0 = [];
sftcAll_90 = [];

for i = 1:length(rst45_135)-1

    sfdum0 = rst0_90{i}.sftc{1};
    sfdum90 = rst0_90{i}.sftc{2};
    sfdum45 = rst45_135{i}.sftc{1};
    sfdum135 = rst45_135{i}.sftc{2};
    
    if size(sfdum45,2) == 5
        sfdum45 = [sfdum45 NaN*ones(size(sfdum45,1),1)];
        sfdum135 = [sfdum135 NaN*ones(size(sfdum135,1),1)];
        sfdum0 = [sfdum0 NaN*ones(size(sfdum0,1),1)];
        sfdum90 = [sfdum90 NaN*ones(size(sfdum90,1),1)];
    end
    
    sftcAll_45 = [sftcAll_45; sfdum45];
    sftcAll_135 = [sftcAll_135; sfdum135];
    sftcAll_0 = [sftcAll_0; sfdum0];
    sftcAll_90 = [sftcAll_90; sfdum90]; 
    
end
    
%%

%vars = {'omag', 'opref' ,'dmag' ,'dpref', 'sfpref', 'sfhco', 'sfBW', 'sfBP','sfCoM'} 

for i = 1:length(vars)
    
    eval(['Q.tc.' vars{i} 'Lum = ' vars{i} 'Lum;'])
    eval(['Q.tc.' vars{i} 'Color = ' vars{i} 'Color;'])
end


Q.sftcAll_45 = sftcAll_45;
Q.sftcAll_135 = sftcAll_135;
Q.sftcAll_0 = sftcAll_0;
Q.sftcAll_90 = sftcAll_90;

Q.Mpk =  max(sftcAll_0');
Q.Spk =  max(sftcAll_90');
Q.Colpk = max(sftcAll_135');
Q.Lumpk = max(sftcAll_45');

Q.M_E =  nanmean(sftcAll_0');
Q.S_E =  nanmean(sftcAll_90');
Q.Col_E = nanmean(sftcAll_135');
Q.Lum_E = nanmean(sftcAll_45');


Q.pSpk = Q.Spk./(Q.Spk+Q.Mpk);
Q.pCpk = Q.Colpk./(Q.Colpk+Q.Lumpk);

Q.pS_E = Q.S_E./(Q.S_E+Q.M_E);
Q.pC_E = Q.Col_E./(Q.Col_E+Q.Lum_E);

Q.pC_Tuning = sftcAll_135./(sftcAll_45+sftcAll_135);

pSlohi = [.2 .8];


%%


idB = find(~isnan(Q.tc.sfprefLum.*Q.tc.sfprefColor.*Q.tc.sfCoMLum.*Q.tc.sfCoMColor.*Q.tc.sfhcoLum.*Q.tc.sfhcoColor.*Q.tc.sfBPLum.*Q.tc.sfBPColor));

analyzeColorSFdist(Q.sftcAll_0(idB,:),Q.sftcAll_45(idB,:),Q.sftcAll_90(idB,:),Q.sftcAll_135(idB,:),pSlohi)



%% 

colorVform(Q.tc,Q.pCpk,Q.pSpk,pSlohi,'Color') %Uses responses to S-M
colorVform(Q.tc,Q.pCpk,Q.pSpk,pSlohi,'Lum') %Uses responses to S+M



%% SF vs. color: Scatter plots comparing color and lum direction of tuning parameters

plotSFparamVColor(Q.tc,Q.pCpk,Q.pSpk,pSlohi)

% 
% ret_percyield = length(find(~isnan(kmap_hor)))/length(kmap_hor)
% overallpercyield = length(find(~isnan(kmap_hor(idB))))/length(kmap_hor)
%% Retinotopy vs. Color 

figure, scatter(kmap_vert,Q.pS_E,'k')
ylabel('%S')
xlabel('Vertical retinotopy')
ylim([0 1])

id = find(~isnan(kmap_vert(:).*Q.pS_E(:)));
[r p] = corrcoef(kmap_vert(id),Q.pS_E(id));
r = r(1,2); p = p(1,2);

title(['p = ' num2str(p)])

%% Orientation comparison

oprefLumD = abs(oridiff(oprefLum*pi/180,0))*180/pi; %projector is rotated 90, so 0 is horizontal
oprefColorD = abs(oridiff(oprefColor*pi/180,0))*180/pi;

figure,
subplot(2,3,1), scatter(oprefLumD,oprefColorD,'.k')
hold on, plot([0 90],[0 90]),xlabel('|oripref-90| (lum)'), ylabel('|oripref-90| (color)')
xlim([0 90]), ylim([0 90]), axis square

subplot(2,3,2), scatter(oprefLum,oprefColor,'.k')
hold on, plot([0 180],[0 180]),xlabel('ori pref (lum)'), ylabel('ori pref (color)')
xlim([0 180]), ylim([0 180]), axis square

subplot(2,3,3), scatter(omagColor,omagLum,'.k')
hold on, plot([0 1],[0 1]),xlabel('ori sel (lum)'), ylabel('ori sel (color)')
axis square



oprefSOpsinD = abs(oridiff(oprefSOpsin*pi/180,0))*180/pi;
oprefMOpsinD = abs(oridiff(oprefMOpsin*pi/180,0))*180/pi;

subplot(2,3,4), scatter(oprefSOpsinD,oprefMOpsinD,'.k')
hold on, plot([0 90],[0 90]),xlabel('|oripref-90| (SOpsin)'), ylabel('|oripref-90| (MOpsin)')
xlim([0 90]), ylim([0 90]), axis square

subplot(2,3,5), scatter(oprefSOpsin,oprefMOpsin,'.k')
hold on, plot([0 180],[0 180]),xlabel('ori pref (SOpsin)'), ylabel('ori pref (MOpsin)')
xlim([0 180]), ylim([0 180]), axis square

subplot(2,3,6), scatter(omagMOpsin,omagSOpsin,'.k')
hold on, plot([0 1],[0 1]),xlabel('ori sel (SOpsin)'), ylabel('ori sel (MOpsin)')
axis square

%% Direction comparison

dprefLumD = abs(angle(exp(1i*(dprefLum-0)*pi/180)))*180/pi; %projector is rotated 90, so 0 is horizontal
dprefColorD = abs(angle(exp(1i*(dprefColor-0)*pi/180)))*180/pi;

figure,
subplot(2,3,1), scatter(dprefLumD,dprefColorD,'.k')
hold on, plot([0 180],[0 180]),xlabel('|dirpref-180| (lum)'), ylabel('|dirpref-180| (color)')
xlim([0 180]), ylim([0 180]), axis square

subplot(2,3,2), scatter(dprefLum,dprefColor,'.k')
hold on, plot([0 360],[0 360]),xlabel('dir pref (lum)'), ylabel('dir pref (color)')
xlim([0 360]), ylim([0 360]), axis square

subplot(2,3,3), scatter(omagColor,omagLum,'.k')
hold on, plot([0 1],[0 1]),xlabel('dir sel (lum)'), ylabel('dir sel (color)')
axis square


dprefSOpsinD = abs(angle(exp(1i*(dprefSOpsin-0)*pi/180)))*180/pi; %projector is rotated 90, so 0 is horizontal
dprefMOpsinD = abs(angle(exp(1i*(dprefMOpsin-0)*pi/180)))*180/pi;

subplot(2,3,4), scatter(dprefSOpsinD,dprefMOpsinD,'.k')
hold on, plot([0 180],[0 180]),xlabel('|dirpref-90| (SOpsin)'), ylabel('|dirpref-90| (MOpsin)')
xlim([0 180]), ylim([0 180]), axis square

subplot(2,3,5), scatter(dprefSOpsin,dprefMOpsin,'.k')
hold on, plot([0 360],[0 360]),xlabel('dir pref (SOpsin)'), ylabel('dir pref (MOpsin)')
xlim([0 360]), ylim([0 360]), axis square

subplot(2,3,6), scatter(omagMOpsin,omagSOpsin,'.k')
hold on, plot([0 1],[0 1]),xlabel('dir sel (SOpsin)'), ylabel('dir sel (MOpsin)')
axis square

%% Sf vs. color selectivity

figure,
subplot(2,2,1), scatter(pColor(idB),sfprefLum(idB)/2+sfprefColor(idB)/2,'.k'), ylabel('sf peak (Lum)'), xlabel('%color')
subplot(2,2,2), scatter(pColor(idB),sfprefColor(idB),'.k'),ylabel('sf peak (color)'), xlabel('%color')
subplot(2,2,3), scatter(pColor(idB),sfCoMLum(idB)/2+sfCoMColor(idB)/2,'.k'), ylabel('sf CoM (Lum)'), xlabel('%color')
subplot(2,2,4), scatter(pColor(idB),sfCoMColor(idB),'.k'),ylabel('sf CoM (color)'), xlabel('%color')


id = find(~isnan(pSAll.*pColor));
[r p] = corrcoef(pSAll(id),pColor(id))
figure,scatter(pSAll,pColor,'.')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
xlabel('%S'), ylabel('%Color')



%% Sf vs. retinotopy

figure,
subplot(2,2,1), scatter(pColor(idB),sfprefLum(idB)/2+sfprefColor(idB)/2,'.k'), ylabel('sf peak (Lum)'), xlabel('%color')
subplot(2,2,2), scatter(pColor(idB),sfprefColor(idB),'.k'),ylabel('sf peak (color)'), xlabel('%color')
subplot(2,2,3), scatter(pColor(idB),sfCoMLum(idB)/2+sfCoMColor(idB)/2,'.k'), ylabel('sf CoM (Lum)'), xlabel('%color')
subplot(2,2,4), scatter(pColor(idB),sfCoMColor(idB),'.k'),ylabel('sf CoM (color)'), xlabel('%color')


id = find(~isnan(pSAll.*pColor));
[r p] = corrcoef(pSAll(id),pColor(id))
figure,scatter(pSAll,pColor,'.')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
xlabel('%S'), ylabel('%Color')


%% sf S vs.M 

idB = find(pSAll>=.0 & pSAll<=1); %comparing opponency vs. summing only makes when there is S AND M contribution
%idB = 1:length(sfprefLum);

sfdom = logspace(log10(.0125),log10(.4),5);

figure,
subplot(4,2,1)
loglog(sfprefSOpsin(idB),sfprefMOpsin(idB),'.k'), hold on, loglog([sfdom(1) .3],[sfdom(1) .3],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf peak loc; S Opsin ')
ylabel('Sf peak loc; M Opsin ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
id = find(~isnan(sfprefSOpsin(idB).*sfprefMOpsin(idB)));
[r p] = corrcoef(log2(sfprefSOpsin(idB(id))),log2(sfprefMOpsin(idB(id))));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(4,2,2)
histogram(log2(sfprefSOpsin(idB)./sfprefMOpsin(idB)),[-4:.25:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfprefSOpsin(idB)./sfprefMOpsin(idB)));
[h p] = ttest(log2(sfprefSOpsin(idB)./sfprefMOpsin(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,3)
loglog(sfCoMSOpsin(idB),sfCoMMOpsin(idB),'.k'), hold on, loglog([sfdom(1) .3],[sfdom(1) .3],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf CoM; S Opsin ')
ylabel('Sf CoM; M Opsin ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
id = find(~isnan(sfCoMSOpsin(idB).*sfCoMMOpsin(idB)));
[r p] = corrcoef(log2(sfCoMSOpsin(idB(id))),log2(sfCoMMOpsin(idB(id))));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(4,2,4)
histogram(log2(sfCoMSOpsin(idB)./sfCoMMOpsin(idB)),[-3:.25:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfCoMSOpsin(idB)./sfCoMMOpsin(idB)));
[hyp p] = ttest(log2(sfCoMSOpsin(idB)./sfCoMMOpsin(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,5)
loglog(sfhcoSOpsin,sfhcoMOpsin,'.k'), hold on, loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf high pass cutoff; S Opsin ')
ylabel('Sf high pass cutoff; M Opsin ')
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
axis square
id = find(~isnan(sfhcoSOpsin(idB).*sfhcoMOpsin(idB)));
[r p] = corrcoef(log2(sfhcoSOpsin(idB(id))),log2(sfhcoMOpsin(idB(id))));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(4,2,6)
histogram(log2(sfhcoSOpsin./sfhcoMOpsin),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfhcoSOpsin./sfhcoMOpsin));
[hyp p] = ttest(log2(sfhcoSOpsin./sfhcoMOpsin));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,7)
scatter(sfBPSOpsin(idB),sfBPMOpsin(idB),'.k'), hold on, plot([0 1],[0 1])
xlabel('Sf bandpass; S Opsin ')
ylabel('Sf bandpass; M Opsin ')
axis square
id = find(~isnan(sfBPSOpsin(idB).*sfBPMOpsin(idB)));
[r p] = corrcoef(sfBPSOpsin(idB(id)),sfBPMOpsin(idB(id)));
title(['r=' num2str(round(r(1,2)*100)/100) '; p=' num2str(p(1,2))])

subplot(4,2,8)
h  = histogram(sfBPSOpsin(idB) - sfBPMOpsin(idB),[-1:.1:1],'FaceColor',[1 1 1]); 
%plot(h.BinEdges(1:end-1),cumsum(h.Values/sum(h.Values)))
set(gca,'TickDir','out')
mu = nanmedian(sfBPSOpsin(idB) - sfBPMOpsin(idB));
[hyp p] = ttest(sfBPSOpsin(idB) - sfBPMOpsin(idB));
title(['mu=' num2str(mu) ' p=' num2str(p)])
xlabel('x-y')




