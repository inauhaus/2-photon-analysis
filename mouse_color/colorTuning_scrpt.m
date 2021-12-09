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

%% Load animal
anim = 'ny6';
expt = 'u003_006';

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')
getCellStats

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories


plot2pColorTuning  %polar plots of color tuning 

%% Load animal
anim = 'rk5';
expt = 'u001_010';

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')
getCellStats

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories


plot2pColorTuning  %polar plots of color tuning 
