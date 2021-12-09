%Script for analyzing kalatsky retinotopy with 2photon data

pF0


global G_handles Analyzer cellS maskS

set(G_handles.epistart,'String','0');  %Frame start in ms (to average)
set(G_handles.epistop,'String','12100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'f:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';



global G_handles

%%
 
clear anim expt maskEx

%Initialize structures with the first animal

% anim{1} = 'nr4';
% expt{1} = 'u006_015';
% maskEx{1} = 'u006_015';

anim{1} = 'ti4';
expt{1} = 'u005_001';
maskEx{1} = 'u005_002';

anim{2} = 'nx8';
expt{2} = 'u001_010';
maskEx{2} = 'u001_010';

anim{3} = 'nx8';
expt{3} = 'u002_005';
maskEx{3} = 'u002_005';

anim{4} = 'nt2';
expt{4} = 'u003_003';
maskEx{4} = 'u003_003';

anim{5} = 'ny6';
expt{5} = 'u002_003';
maskEx{5} = 'u002_003';



eid = 1


%load expt
set(G_handles.datadir,'string',[dataRoot anim{eid} '\' anim{eid} '_' expt{eid}(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim{eid} '\' anim{eid} '_' expt{eid}(1:8)])
Gsetdirectories

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim{eid} '_' maskEx{eid}(1:8)];
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim{eid} '_' expt{eid}(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

%f1 = f1meanimage;  %Build F1 images (takes the longest)
global Analyzer
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


%% Initialize structures with the first animal

global cellS maskS Analyzer

%Nice anatomy with lots of cells in mask, but most seem poorly tuned.
%Result of color vs. sf is consistent

anim = 'rl0';

expt = 'u002_006';
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_001'; %Kalatsky
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_6'];
load(tracepath,'cellS')
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

%f1 = f1meanimage;  %Build F1 images (takes the longest)
%First set dirs and hit 'Compute F0 images' within pF0
f1 = CondF1_cellMask2(.5);

[kmap_hor kmap_vert] = processkret_cellMask(f1);  %Make maps to plot, delete L if no smoothing

projectorAdjustment = 1;
if projectorAdjustment    
    kmap_horx = kmap_hor;
    kmap_vertx = kmap_vert;
    kmap_hor = kmap_vertx;
    kmap_vert = kmap_horx; %This should be negative if run with LCD rotated 90 clockwise
end

plotMaskedRetinotopy(kmap_hor,kmap_vert) 
