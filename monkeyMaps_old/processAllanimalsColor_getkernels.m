pRev

global G_RChandles G_handles

set(G_RChandles.kernelLength,'string','[-500 2000]');
set(G_RChandles.LPflag,'value',0);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);
set(G_RChandles.blankNorm,'value',0);

set(G_handles.datadir,'string','C:\2p_data\')
set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

set(G_handles.fastMotionFlag,'value',1)
set(G_handles.slowMotionFlag,'value',1)
set(G_handles.searchRange,'string','5')
set(G_handles.blockSize,'string','30')

set(G_handles.F0flag,'value',0)
set(G_handles.cellmaskflag,'value',1)
set(G_handles.loadSyncs,'value',1)

% set(G_handles.datadir,'string','g:\2p_data\')
% set(G_handles.analyzedir,'string','g:\AnalyzerFiles\')

%% ab2; u000_011

global maskS cellS

anim = 'ab2';
expt = 'u000_011';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 

%% ab2; u000_057

global maskS cellS

anim = 'ab2';
expt = 'u000_057';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 


%% ab2; u000_089

global maskS cellS

anim = 'ab2';
expt = 'u000_089';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 

%% ab2; u002_018

global maskS cellS

anim = 'ab2';
expt = 'u002_018';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 

%% ab3; u002_013

global maskS cellS

anim = 'ab3';
expt = 'u002_013';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 


%% ab3; u002_042

global maskS cellS

anim = 'ab3';
expt = 'u002_042';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS') 

%Get traces
%processButton

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);

%save traces and kernels
save(tracepath,'cellS') 


% %% aa9; u001_026
% 
% global maskS cellS
% 
% anim = 'aa9';
% expt = 'u001_026';
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% 
% maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
% maskpath = [maskroot anim '_' expt];
% load(maskpath,'maskS') 
% 
% traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
% tracepath = [traceroot anim '_' expt '_cellS'];
% load(tracepath,'cellS') 
% 
% %Get traces
% %processButton
% 
% %make Kernels
% set(G_RChandles.dropTrials,'string','[]')
% hh = makeTemporalfilter;
% trialdom = 1:1:getnotrials;
% eval(['dT = ' get(G_RChandles.dropTrials,'string')])
% trialdom(dT) = [];
% Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);
% 
% %save traces and kernels
% save(tracepath,'cellS') 
% 
% set(G_handles.blockSize,'string','30') %return to default
% 
