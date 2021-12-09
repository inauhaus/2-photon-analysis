pF0
pRev
%%
global TimingInfo G_handles cellS maskS

anim = 'fry';
expt = 'u008_001';

expt = 'u006_006';

%expt = 'u008_003';

%expt = 'u005_005'; %randori

expt = 'u005_004'; %cartoons

datadir = ['f:\2p_data\' anim '\' anim '_' expt(2:end)];
anadir = ['c:\AnalyzerFiles\' anim '\' anim '_' expt];

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS')

traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

load(['C:\Slave_Files\' anim '\' anim '_' expt(2:end)])


SaccadeTriggeredResponse