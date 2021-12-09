pF0

%initialize the Gui

global G_handles Analyzer cellS maskS

set(G_handles.epistart,'String','100');  %Frame start in ms (to average)
set(G_handles.epistop,'String','1100'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-50');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'e:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';

%% Load raindropper expt
global cellS maskS

anim = 'ny6';

expt = 'u003_006'; %use mask from other expt
%expt = 'u003_009'; %use mask from other expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u003_009'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 3_6'];
%tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);


%% Load raindropper expt
global cellS maskS

anim = 'ny6';

expt = 'u002_006'; %use mask from other (or same) expt
%expt = 'u002_008'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_006'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_8'];
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'nz6';

expt = 'u002_013'; %use mask from other (or same) expt
%expt = 'u002_011'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_011'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_13'];
%tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);


%% Load raindropper expt
global cellS maskS

anim = 'nz6';

expt = 'u001_007'; %use mask from other (or same) expt
%expt = 'u001_005'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_005'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_7'];
%tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rc3';

expt = 'u005_005'; %use mask from other (or same) expt
expt = 'u005_004'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u005_004'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 5_5'];
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rk5';

expt = 'u001_010'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_008'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_10'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rl0';

expt = 'u001_005'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_009'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_5'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rk3';

expt = 'u002_009'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_013'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_9'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);


%% Load raindropper expt
global cellS maskS

anim = 'rk3';

expt = 'u001_010'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_011'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_10'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rk4';

expt = 'u001_003'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u001_004'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 1_3'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rk5';

expt = 'u002_002'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_010'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS_aligned to 2_2'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);

%% Load raindropper expt
global cellS maskS

anim = 'rl0';

expt = 'u002_012'; %use mask from other (or same) expt
maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
load(maskpath,'maskS') 

expt = 'u002_012'; %raindropper
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt '_cellS'];
load(tracepath,'cellS')

set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

getCellStats

xsize = getparam('x_size');
posdom = getdomain('s_phase')*xsize/360;
KalLocs = 30;  % RF locations from Kalatsky
anatCorrection = 30;
posdom = posdom+KalLocs-anatCorrection;

rst = episodicRaindropper(posdom);



