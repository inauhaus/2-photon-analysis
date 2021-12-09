pF0

%initialize the Gui

global G_handles Analyzer cellS maskS

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%% Load raindropper expt
global cellS maskS

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','1200'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-100');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

anim = 'ny6';
Mexpt = 'u003_006'; %mask expt
Dexpt = 'u003_009'; %data expt

% anim = 'ny6';
% Mexpt = 'u002_006'; %mask expt
% Dexpt = 'u002_006'; %data expt

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' Mexpt(1:8)];
load(maskpath,'maskS') 

traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' Dexpt '_cellS_aligned to 3_6'];
%tracepath = [traceroot anim '_' Dexpt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' Dexpt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' Dexpt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

figure,scatter(2*rst.EnvM.param(:,2),2*rst.EnvS.param(:,2)), 
hold on, 
plot([0 40],[0 40]),
xlabel('RF width (2sigma) Mopsin'), ylabel('RF width (2sigma) Sopsin'), axis square
xlim([0 40]), ylim([0 40])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Drifting color expt

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-400');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

anim = 'ny6';
expt = 'u003_006';
%expt = 'u002_004';

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories
getCellStats

X = plot2pColorTuning;


%%
%Compute parameters from grating experiment
pColor =  X.pC;
pS =  X.pS;
Cdir = X.dir;

%Compute parameters from rain experiment
ES = rst.EnvS.param(:,3);
EM = rst.EnvM.param(:,3);
% ES = sum(rst.EnvS.ffit,2);
% EM = sum(rst.EnvM.ffit,2);
pSE = ES./(EM+ES);

Son = sum(rst.sON.ffit,2);
Soff = sum(rst.sOFF.ffit,2);
Mon = sum(rst.mON.ffit,2);
Moff = sum(rst.mOFF.ffit,2);

Son = rst.sON.param(:,3);
Soff = rst.sOFF.param(:,3);
Mon = rst.mON.param(:,3);
Moff = rst.mOFF.param(:,3);

colorVec = Mon-Moff + 1i*(Son-Soff);
CdirE = angle(colorVec')*180/pi;
CdirE(find(CdirE<0)) =  CdirE(find(CdirE<0)) + 180;

id = find(isnan(Cdir.*CdirE));
Cdir(id) = [];
CdirE(id) = [];

figure,scatter(Cdir,CdirE), 
figure,scatter(pS,pSE), xlabel('%S drifting gratings'), ylabel('%S rain')
hold on, plot([0 1],[0 1])

figure,hist(X.pC,[0:.1:1]), xlabel('color/(color+lum)')

[theta] = circCorr(Cdir,CdirE,180)
