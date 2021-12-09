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

%% ab2; u000_017

global maskS cellS

anim = 'ab2';
expt = 'u000_017';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 

%% ab2; u000_070

global maskS cellS

anim = 'ab2';
expt = 'u000_070';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 


%% ab2; u001_006

global maskS cellS

anim = 'ab2';
expt = 'u001_006';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[13:40]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 

%% aa9; u001_031

global maskS cellS

anim = 'aa9';
expt = 'u001_031';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 

%% ab3; u002_017

global maskS cellS

anim = 'ab3';
expt = 'u002_017';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 


%% ab8; u001_035

global maskS cellS

anim = 'ab8';
expt = 'u001_035';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 


%% ab4; u000_027

set(G_handles.blockSize,'string','60') %very little heart beat motion, so make this a bit longer 

global maskS cellS

anim = 'ab4';
expt = 'u000_027';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[21:30]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 

set(G_handles.blockSize,'string','30') %return to default

%% ab3; u002_028

global maskS cellS

anim = 'ab3';
expt = 'u002_028';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
load(maskpath,'maskS') 

traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];
%load(tracepath,'cellS') 

%Get traces
try
    processButton
end

%make Kernels
set(G_RChandles.dropTrials,'string','[]')
hh = makeTemporalfilter;
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
Ggetrandposkernel2(cellS.cellMat,trialdom,hh)

%save traces and kernels
save(tracepath,'cellS') 

