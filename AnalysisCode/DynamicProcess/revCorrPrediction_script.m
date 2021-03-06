
pRev

global G_RChandles G_handles

set(G_RChandles.kernelLength,'string','[-500 2000]');
set(G_RChandles.LPflag,'value',1);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);
set(G_RChandles.blankNorm,'value',0);

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

set(G_handles.datadir,'string','e:\2p_data\')
set(G_handles.analyzedir,'string','e:\2p_data\AnalyzerFiles\')

%%

global maskS cellS idExamp

anim = 'ab2';
expt = 'u000_014';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
%traceroot = 'e:\Beta\cellTraces\monkey\new2\';
%traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS')

%idExamp = [6 27]
%idExamp = [79 59 33]
%idExamp = [6 34 54 79];  %ab3 002_041, ab2  000_068
%idExamp = [9 54 78 83];  %aa9 001_027
idExamp = [3 31 54];  %ab2 000_014


%% mouse Gcamp

global maskS cellS idExamp

anim = 'zc1';
expt = 'u000_025';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\mouse\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\2ph_code\Beta\cellTraces\mouse\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS')

%idExamp = [6 27]
%idExamp = [79 59 33]
%idExamp = [6 34 54 79];  %ab3 002_041, ab2  000_068
%idExamp = [9 54 78 83];  %aa9 001_027
idExamp = [3 31 54];  %ab2 000_014

%%

kern = revCorrPrediction(cellS.cellMat,1:30)

