pF0

global G_handles
global maskS cellS

set(G_handles.epistart,'String','0');  %Frame start in ms (to average)
set(G_handles.epistop,'String','2000'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-500');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'value',1); %Frame stop in ms (to average)

dataRoot = 'g:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%%

anim = 'np1';
expt = 'u001_004';

%load expt
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

%Process Data
[kernPop{1} kernSigPop{1}] = getTCmat_greenblue;
[percS{1} rfit{1}] = getOpsinInput(kernPop{1});

plotOpsinMap(percS{1},rfit{1})

%%%%



