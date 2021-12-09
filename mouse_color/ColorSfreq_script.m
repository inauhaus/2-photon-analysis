pF0

%initialize the Gui

global G_handles Analyzer cellS maskS

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

clear rst0_90 rst45_135
%% Load color x sf x ori expt
global cellS maskS idExamp

anim = 'rk5'; 

eid = 1;

idExamp = [3 21 75 81 92];
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_010'; 
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_010'; %S vs M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%
%Use mask from 001_010
expt = 'u001_011'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_10'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);

%% Load color x sf x ori expt
%weak response to color direction, high %S, coincidence?
global cellS maskS

anim = 'rk7';
eid = 2;

idExamp = [3 21 ];
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_006';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_006'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u001_007'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_6'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);

%% Load color x sf x ori expt 
global cellS maskS
idExamp = [19 21 25 49 58 78 98];
%good %S map
%nice anatomy. lots of cells.
%Very strong Bandpass vs color effect
%oripref has little correspondence

anim = 'rl0';
eid = 3;
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_005';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_005'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u001_006'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_5'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);

plotBothColorDirs(rst0_90{eid},rst45_135{eid})

%% Load color x sf x ori expt
global cellS maskS

anim = 'rl1';
eid = 4;
idExamp = [];
%not many cells. Kinda crappy.
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u002_006';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_006'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u002_007'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_6'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

% rst45_135{eid} = ColorSfreq2(1,2);

%% Load color x sf x ori expt
global cellS maskS

anim = 'rk3'; %good
eid = 5;
%good %S map and SF vs color effect
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u002_009';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_010'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_9'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u002_011'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_9'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);

%% Load color x sf x ori expt
global cellS maskS
%perfect %s map, 
%sfpeak bias is in the other direction 
%sfCoM and BP are consistent
anim = 'rl2';
eid = 6;
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u000_005';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u000_006'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 0_5'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u000_007'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 0_5'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);

%% Load color x sf x ori expt
global cellS maskS

anim = 'rk4'; %This experiment has all 4 color direction together
eid = 7;
%ok; sf bias is in correct direction.  Consistently low-pass in color
%direction
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_003';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_005'; %S,S+M,M,S-M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_3'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

c1 = 1; c2 = 3;
rst0_90{eid} = ColorSfreq2(c1,c2);

c1 = 2; c2 = 4;
rst45_135{eid} = ColorSfreq2(c1,c2);

%% Load color x sf x ori expt
global cellS maskS

anim = 'rk5';
eid = 8;
%Not many good cells. Mostly sf vs. color seems mostly separable

%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u002_002';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_004'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_2'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u002_005'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_2'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);


%% Load color x sf x ori expt
global cellS maskS

anim = 'rl0';
eid = 9;
%Nice anatomy with lots of cells in mask, but most seem poorly tuned.
%Result of color vs. sf is consistent
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u002_006';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_004'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_6'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst0_90{eid} = ColorSfreq2(1,2);

%%%%%%Now get 45 and 135%%%%%%%%

expt = 'u002_006'; %color vs lum
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

rst45_135{eid} = ColorSfreq2(1,2);



%% Population measures
sfprefColor = [];
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
pS = [];
pColor = [];

for i = 1:length(rst45_135)
%for i = 1:3
    
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
    end
    
end

idB = find(pS>.0 & pS<1); %comparing opponency vs. summing only makes when there is S AND M contribution

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

% subplot(4,2,5)
% loglog(sfhcoLum,sfhcoColor,'.k'), hold on, loglog([sfdom(1) sfdom(end) ],[sfdom(1) sfdom(end)],'r')
% set(gca,'XTick', sfdom(1:2:end), 'YTick', sfdom(1:2:end))
% xlabel('Sf high pass cutoff; color 1 ')
% ylabel('Sf high pass cutoff; color 2 ')
% xlim([sfdom(1) sfdom(end)]), ylim([sfdom(1) sfdom(end)])
% axis square
% 
% subplot(4,2,6)
% histogram(log2(sfhcoLum./sfhcoColor),[-4:.5:4],'FaceColor',[1 1 1])
% set(gca,'TickDir','out')
% mu = nanmean(log2(sfhcoLum./sfhcoColor));
% [hyp p] = ttest(log2(sfhcoLum./sfhcoColor));
% title(['geomu=' num2str(2.^mu) ' p=' num2str(p)])
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

