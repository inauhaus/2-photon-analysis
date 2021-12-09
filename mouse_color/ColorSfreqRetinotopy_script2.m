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

clear rst0_90 rst45_135

%%
animAll{1} = 'rk5';
expt_mask{1} = 'u001_010';
expt_0v90{1} = 'u001_010';
expt_45v135{1} = 'u001_011';
expt_Kal{1} = 'u001_002';
idExampAll{1} = [3 21 75 81 92];

animAll{2} = 'rk7';
expt_mask{2} = 'u001_006';
expt_0v90{2} = 'u001_006';
expt_45v135{2} = 'u001_007';
expt_Kal{2} = 'u001_000';
idExampAll{2} = [3 21];

animAll{3} = 'rl0';
expt_mask{3} = 'u001_005';
expt_0v90{3} = 'u001_005';
expt_45v135{3} = 'u001_006';
expt_Kal{3} = 'u001_001';
idExampAll{3} = [19 21 25 49 58 78 98];

animAll{4} = 'rl1';
expt_mask{4} = 'u002_006';
expt_0v90{4} = 'u002_006';
expt_45v135{4} = 'u002_007';
expt_Kal{4} = 'u002_003';
idExampAll{4} = [];

animAll{5} = 'rk3';
expt_mask{5} = 'u002_009';
expt_0v90{5} = 'u002_010';
expt_45v135{5} = 'u002_011';
expt_Kal{5} = 'u002_004';
idExampAll{5} = [];

animAll{6} = 'rl2';
expt_mask{6} = 'u000_005';
expt_0v90{6} = 'u000_006';
expt_45v135{6} = 'u000_007';
expt_Kal{6} = 'u000_001';
idExampAll{6} = [];

animAll{7} = 'rk5';
expt_mask{7} = 'u002_002';
expt_0v90{7} = 'u002_004';
expt_45v135{7} = 'u002_005';
expt_Kal{7} = 'u002_001';
idExampAll{7} = [];

animAll{8} = 'rl0';
expt_mask{8} = 'u002_006';
expt_0v90{8} = 'u002_004';
expt_45v135{8} = 'u002_006';
expt_Kal{8} = 'u002_001';
idExampAll{8} = [];

animAll{9} = 'rk4';
expt_mask{9} = 'u001_003'; 
expt_0v90{9} = 'u001_005';
expt_45v135{9} = 'u001_002';
expt_Kal{9} = 'u001_002';
idExampAll{9} = [];
%% Load color x sf x ori expt

global cellS maskS idExamp

exdom = 1:8;
%exdom = 4;
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
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    
    try
        Gsetdirectories %Load raw experiment and analyzer file
    catch
        loadAnalyzer %Just load analyzer file
    end
    %%%%%%%%%%%
    
    getCellStats
    
    rst0_90{ex} = ColorSfreq2(1,2);
    
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
    
    rst45_135{ex} = ColorSfreq2(1,2);
    
    plotBothColorDirs(rst0_90{ex},rst45_135{ex})
    
    
    %%%%%%%%%%%%%%%%%Now get Kalatsky info and compare to S/M %%%%%%%%%%%%%%%%
    expt = expt_Kal{ex}; %Kalatsky
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS_aligned to ' expt_mask{ex}(2:end)];
    load(tracepath,'cellS')
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    
    getCellStats
    
    %f1 = f1meanimage;  %Build F1 images (takes the longest)
    %First set dirs and hit 'Compute F0 images' within pF0
    f1 = CondF1_cellMask2(.3);
    
    [kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing
    
    projectorAdjustment = 1;
    
    if projectorAdjustment
        kmap_horx = kmap_hor;
        kmap_vertx = kmap_vert;
        kmap_hor = kmap_vertx;
        kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
    end
    
    plotMaskedRetinotopy(kmap_hor,kmap_vert)
    
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

%%


%% Load color x sf x ori expt
global cellS maskS

anim = 'rk4'; %This experiment has all 4 color direction together
eid = 9;
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
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

%f1 = f1meanimage;  %Build F1 images (takes the longest)
%First set dirs and hit 'Compute F0 images' within pF0
f1 = CondF1_cellMask2(.3);

[kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing

projectorAdjustment = 1;

if projectorAdjustment    
    kmap_horx = kmap_hor;
    kmap_vertx = kmap_vert;
    kmap_hor = kmap_vertx;
    kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
end

plotMaskedRetinotopy(kmap_hor,kmap_vert) 

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
rst45_135{i}


sfpref_135 = [];
sfpref_0 = [];
sfpref_90 = [];
sfpref_135 = [];

sfprefMopsin = [];
sfprefLum = [];
sfprefCoMColor = [];
sfhcoLum = [];
sfhcoColor = [];
sfprefCoMLum = [];
sfBPColor = [];
sfBPLum = [];
oprefL = [];
oprefC = [];
omagL = [];
omagC = [];
kmap_vert = [];
kmap_hor = [];
pS = [];
pColor = [];

%for i = 1:length(rst45_135)
for i = 1:8
    
    if ~isempty(rst45_135{i})
        sfprefColor = [sfprefColor rst45_135{i}.sfpref{2}];
        sfprefLum = [sfprefLum rst45_135{i}.sfpref{1}];
        
        sfprefCoMColor = [sfprefCoMColor rst45_135{i}.sfCoM{2}];
        sfprefCoMLum = [sfprefCoMLum rst45_135{i}.sfCoM{1}];
        
        sfhcoColor = [sfhcoColor rst45_135{i}.sfhco{2}];
        sfhcoLum = [sfhcoLum rst45_135{i}.sfhco{1}];
        
        sfBPColor = [sfBPColor rst45_135{i}.sfBP{2}];
        sfBPLum = [sfBPLum rst45_135{i}.sfBP{1}];
        
        oprefC = [oprefC rst45_135{i}.opref{2}];
        oprefL = [oprefL rst45_135{i}.opref{1}];
        
        omagC = [omagC rst45_135{i}.omag{2}];
        omagL = [omagL rst45_135{i}.omag{1}];
        
        pColor = [pColor rst45_135{i}.pS];
        pS = [pS rst0_90{i}.pS];
        
        kmap_hor = [kmap_hor Kret{i}.kmap_hor'];
        kmap_vert = [kmap_vert Kret{i}.kmap_vert'];
        
    end
    
end


idB = find(pS>.0 & pS<1); %comparing opponency vs. summing only makes when there is S AND M contribution

%%

sfdom = logspace(log10(.0125),log10(.4),5);

figure,
subplot(4,2,1)
loglog(sfprefLum(idB),sfprefColor(idB),'.k'), hold on, loglog([sfdom(1) .3],[sfdom(1) .3],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf peak loc; color 1 ')
ylabel('Sf peak loc; color 2 ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])

subplot(4,2,2)
histogram(log2(sfprefLum(idB)./sfprefColor(idB)),[-4:.25:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfprefLum(idB)./sfprefColor(idB)));
[h p] = ttest(log2(sfprefLum(idB)./sfprefColor(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,3)
loglog(sfprefCoMLum(idB),sfprefCoMColor(idB),'.k'), hold on, loglog([sfdom(1) .3],[sfdom(1) .3],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf CoM; color 1 ')
ylabel('Sf CoM; color 2 ')
axis square
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])

subplot(4,2,4)
histogram(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)),[-3:.25:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)));
[hyp p] = ttest(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,5)
loglog(sfhcoLum,sfhcoColor,'.k'), hold on, loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'r')
set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
xlabel('Sf high pass cutoff; color 1 ')
ylabel('Sf high pass cutoff; color 2 ')
xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
axis square

subplot(4,2,6)
histogram(log2(sfhcoLum./sfhcoColor),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfhcoLum./sfhcoColor));
[hyp p] = ttest(log2(sfhcoLum./sfhcoColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
% xlabel('log(x)-log(y)')

subplot(4,2,7)
scatter(sfBPLum(idB),sfBPColor(idB),'.k'), hold on, plot([0 1],[0 1])
xlabel('Sf bandpass; color 1 ')
ylabel('Sf bandpass; color 2 ')
axis square

subplot(4,2,8)
h  = histogram(sfBPLum(idB) - sfBPColor(idB),[-1:.1:1],'FaceColor',[1 1 1]); 
%plot(h.BinEdges(1:end-1),cumsum(h.Values/sum(h.Values)))
set(gca,'TickDir','out')
mu = nanmedian(sfBPLum(idB) - sfBPColor(idB));
[hyp p] = ttest(sfBPLum(idB) - sfBPColor(idB));
title(['mu=' num2str(mu) ' p=' num2str(p)])
xlabel('x-y')

%% SF vs. retinotopy

sfdom = logspace(log10(.0125),log10(.4),5);

sfprefmu = (sfprefLum .* sfprefColor .* sfprefLum .* sfprefColor).^(1/4);

figure,
subplot(4,2,1)
semilogx(sfprefmu(idB),kmap_vert(idB),'.k')
axis square
set(gca,'XTick', sfdom(1:2:end), 'YTick', [-45 0 45])
xlabel('Sf peak loc')
ylabel('retinotopy')
xlim([sfdom(1) sfdom(end)]), 
ylim([-45 45])

subplot(4,2,2)
histogram(log2(sfprefLum(idB)./sfprefColor(idB)),[-4:.25:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfprefLum(idB)./sfprefColor(idB)));
[h p] = ttest(log2(sfprefLum(idB)./sfprefColor(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,3)
semilogx(sfprefCoMLum(idB),kmap_vert(idB),'.k')
axis square
set(gca,'XTick', sfdom(1:2:end), 'YTick', [-45 0 45])
xlabel('Sf CoM')
ylabel('retinotopy')
xlim([sfdom(1) sfdom(end)]), 
ylim([-45 45])

subplot(4,2,4)
histogram(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)),[-3:.25:3],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)));
[hyp p] = ttest(log2(sfprefCoMLum(idB)./sfprefCoMColor(idB)));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
xlabel('log(x)-log(y)')

subplot(4,2,5)
semilogx(sfhcoLum(idB),kmap_vert(idB),'.k')
axis square
set(gca,'XTick', sfdom(1:2:end), 'YTick', [-45 0 45])
xlabel('Sf hp cutoff')
ylabel('retinotopy')
xlim([sfdom(1) sfdom(end)]), 
ylim([-45 45])

subplot(4,2,6)
histogram(log2(sfhcoLum./sfhcoColor),[-4:.5:4],'FaceColor',[1 1 1])
set(gca,'TickDir','out')
mu = nanmean(log2(sfhcoLum./sfhcoColor));
[hyp p] = ttest(log2(sfhcoLum./sfhcoColor));
title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
% xlabel('log(x)-log(y)')

subplot(4,2,7)
semilogx(sfBPLum(idB),kmap_vert(idB),'.k')
axis square
set(gca,'XTick', sfdom(1:2:end), 'YTick', [-45 0 45])
xlabel('Sf BP factor')
ylabel('retinotopy')
xlim([sfdom(1) sfdom(end)]), 
ylim([-45 45])

subplot(4,2,8)
h  = histogram(sfBPLum(idB) - sfBPColor(idB),[-1:.1:1],'FaceColor',[1 1 1]); 
%plot(h.BinEdges(1:end-1),cumsum(h.Values/sum(h.Values)))
set(gca,'TickDir','out')
mu = nanmedian(sfBPLum(idB) - sfBPColor(idB));
[hyp p] = ttest(sfBPLum(idB) - sfBPColor(idB));
title(['mu=' num2str(mu) ' p=' num2str(p)])
xlabel('x-y')


%% Retinotopy vs. Color 

figure, scatter(kmap_vert,pS,'k')
ylabel('%S')
xlabel('Vertical retinotopy')
ylim([0 1])

id = find(~isnan(kmap_vert(:).*pS(:)));
[r p] = corrcoef(kmap_vert(id),pS(id));
r = r(1,2); p = p(1,2);

title(['p = ' num2str(p)])

%% Orientation comparison
% 
% id = find(oprefL>90);
% oprefL(id) = oprefL(id)-90;
% 
% id = find(oprefC>90);
% oprefC(id) = oprefC(id)-90;

edg = 0:10:180;
h = histogram(oprefL,edg);
hL = h.Values; close;
h = histogram(oprefC,edg);
hC = h.Values; close;
hdom = edg(1:end-1)-(edg(1)-edg(1))

figure, 
subplot(2,2,1)
histogram(oprefL,edg,'FaceColor',[1 1 1]); 
title('Luminance condition')
subplot(2,2,3)
histogram(oprefC,edg,'FaceColor',[1 1 1]); 
title('Color condition')
xlabel('oripref')

subplot(2,2,2)
stairs(cumsum(hL))
hold on
stairs(cumsum(hC))


figure,
subplot(1,2,1), scatter(oprefL,oprefC,'.k')
hold on, plot([0 180],[0 180]),xlabel('ori pref (lum)'), ylabel('ori pref (color)')
xlim([0 180]), ylim([0 180]), axis square
subplot(1,2,2), scatter(omagC,omagL,'.k')
hold on, plot([0 1],[0 1]),xlabel('ori sel (lum)'), ylabel('ori sel (color)')
axis square

%% Sf vs. color selectivity

figure,
subplot(2,2,1), scatter(pColor(idB),sfprefLum(idB)/2+sfprefColor(idB)/2,'.k'), ylabel('sf peak (Lum)'), xlabel('%color')
subplot(2,2,2), scatter(pColor(idB),sfprefColor(idB),'.k'),ylabel('sf peak (color)'), xlabel('%color')
subplot(2,2,3), scatter(pColor(idB),sfprefCoMLum(idB)/2+sfprefCoMColor(idB)/2,'.k'), ylabel('sf CoM (Lum)'), xlabel('%color')
subplot(2,2,4), scatter(pColor(idB),sfprefCoMColor(idB),'.k'),ylabel('sf CoM (color)'), xlabel('%color')


id = find(~isnan(pS.*pColor));
[r p] = corrcoef(pS(id),pColor(id))
figure,scatter(pS,pColor,'.')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
xlabel('%S'), ylabel('%Color')



%% Sf vs. retinotopy

figure,
subplot(2,2,1), scatter(pColor(idB),sfprefLum(idB)/2+sfprefColor(idB)/2,'.k'), ylabel('sf peak (Lum)'), xlabel('%color')
subplot(2,2,2), scatter(pColor(idB),sfprefColor(idB),'.k'),ylabel('sf peak (color)'), xlabel('%color')
subplot(2,2,3), scatter(pColor(idB),sfprefCoMLum(idB)/2+sfprefCoMColor(idB)/2,'.k'), ylabel('sf CoM (Lum)'), xlabel('%color')
subplot(2,2,4), scatter(pColor(idB),sfprefCoMColor(idB),'.k'),ylabel('sf CoM (color)'), xlabel('%color')


id = find(~isnan(pS.*pColor));
[r p] = corrcoef(pS(id),pColor(id))
figure,scatter(pS,pColor,'.')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
xlabel('%S'), ylabel('%Color')
