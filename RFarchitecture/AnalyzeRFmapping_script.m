%% Open the GUI and set some timing parameters

pF0

%Set the paths for this new code
path('C:\2pScanboxAnalysis\RFarchitecture',path)

global G_handles Analyzer cellS maskS idExamp

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-500');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'd:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%% ID the experiments to analyze

% anim = 'xl2';
% exptRet = 'u017_002'; %Episodic retinotopy experiment
% exptSF = 'u017_009';  %Spatial frequency tuning experiment
% exptTemp = 'u017_005';  %This is the experiment used to create the mask

% anim = 'xl2';
% exptRet = 'u018_003'; %Episodic retinotopy experiment
% exptSF = 'u018_006';  %Spatial frequency tuning experiment
% exptTemp = 'u018_005';  %This is the experiment used to create the mask

% anim = 'xl2';
% exptRet = 'u020_004'; %Episodic retinotopy experiment
% exptSF = 'u018_006';  %Spatial frequency tuning experiment
% exptTemp = 'u018_005';  %This is the experiment used to create the mask

% anim = 'ws5';
% exptRet = 'u013_002'; %Episodic retinotopy experiment
% exptSF = 'u013_005';  %Spatial frequency tuning experiment
% exptTemp = 'u013_003';  %This is the experiment used to create the mask
% 
% anim = 'xb7';
% exptRet = 'u014_002'; %Episodic retinotopy experiment
% exptSF = 'u014_006';  %Spatial frequency tuning experiment
% exptTemp = 'u014_001';  %This is the experiment used to create the mask


%% Load the RF mapping experiment
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' exptRet(2:end) ''])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' exptRet(1:8)])
Gsetdirectories

%%Load mask of the template experiment%%%%%%%
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' exptTemp(2:8)];
load(maskpath,'maskS')

%%Load traces for episodic retinotopy experiment%%%%%%%
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' exptRet(1:8) '_cellS_alignedto' num2str(exptTemp(3:4)) '_' num2str(exptTemp(8))];
load(tracepath,'cellS')


%%%%%%%Compute RF parameters (fit curves, etc)%%%%%%%%%

RF = getRFsizepos;  %RF is a strucutre with the tuning fits and parameters

CoMyx = getCellPositions;  % [y x] location of each cell

%%%%%%%Plot RF position and size maps%%%%%%%

vThresh = 0.7;  %Threshold for fit quality.  Variance accounted for.

plotRFpositions(CoMyx,RF,vThresh)

plotRFsizes(CoMyx,RF,vThresh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SF tuning analysis

%%Load traces%%%%%%%
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' exptSF(1:8) '_cellS_alignedto' num2str(exptTemp(3:4)) '_' num2str(exptTemp(8))];
load(tracepath,'cellS')

SF = getSFtuning;  

vThresh = 0.5;
plotSF(CoMyx,SF,vThresh)

