pF0

global G_handles Analyzer cellS maskS

dataRoot = 'h:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

%%

anim = 'rl0';
%Nice anatomy with lots of cells in mask, but most seem poorly tuned.
%Result of color vs. sf is consistent
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
rst = ColorSfreq2(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expt = 'u001_001'; %Kalatsky
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 001_005'];
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

figure, scatter(kmap_vert,rst.pS,'k')
ylabel('%S')
xlabel('Vertical retinotopy')
ylim([0 1])

id = find(~isnan(kmap_vert(:).*rst.pS(:)));
[r p] = corrcoef(kmap_vert(id),rst.pS(id));
r = r(1,2); p = p(1,2);

title(['p = ' num2str(p)])

%%
%%

anim = 'rk5';
idExamp = [3 21 75 81 92];
%Nice anatomy with lots of cells in mask, but most seem poorly tuned.
%Result of color vs. sf is consistent
%%%%%%First get 0 and 90%%%%%%%%%%%%

expt = 'u001_010';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_010'; %S vs. M
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats
rst = ColorSfreq2(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

expt = 'u001_002'; %Kalatsky
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 001_010'];
load(tracepath,'cellS')
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

%f1 = f1meanimage;  %Build F1 images (takes the longest)
%First set dirs and hit 'Compute F0 images' within pF0
f1 = CondF1_cellMask2(.4);

[kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing

projectorAdjustment = 1;

if projectorAdjustment    
    kmap_horx = kmap_hor;
    kmap_vertx = kmap_vert;
    kmap_hor = kmap_vertx;
    kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
end

plotMaskedRetinotopy(kmap_hor,kmap_vert) 

figure, scatter(rst.pS,kmap_vert)
