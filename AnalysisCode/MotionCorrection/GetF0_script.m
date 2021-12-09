pF0

%initialize the Gui

global G_handles Analyzer cellS maskS

set(G_handles.epistart,'String','500');  %Frame start in ms (to average)
set(G_handles.epistop,'String','3250'); %Frame stop in ms (to average)
set(G_handles.bstart,'String','-1000');  %Frame start in ms (to average)
set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

dataRoot = 'f:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';


%%
clear animAll expt_mask

id = 1;

%%%%%%%%%%%%%%%
animAll{id} = 'fry';
exptAll{id} = 'u003_001'; %ori8

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u003_002'; %ori2 sf = [.2 .4 .8];  decent responses

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u003_003'; %ori8; sf = [.5 1 2 4]

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u003_009';  %ori [0 45 90]; contrast4

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u003_010'; %ori [0 45 90 135]; contrast 5

%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u004_001'; %ori4;sf2

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u004_002'; %contrast

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u004_003';  %flashed gratings
 
%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u005_004'; 

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u005_005'; 

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u005_006';

%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_006'; %flashed grater

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_007'; %eye pos calibrator

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_008'; %eye pos calibrator

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_009'; %texture vs. noise

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_010'; %Spiderman

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u006_012'; %contrast reversing gratings

%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u007_002'; %ori 8

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u007_005'; %face calibrator

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u007_006'; %face calibrator

%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_001'; %Spiderman

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_002'; %texture v noise

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_003'; %random gratings

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_004'; %periodic retinotopy

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_005'; %contrast reversing gratings
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_006'; %texture v noise

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u008_008'; %random gratings

%%
clear animAll exptAll
animAll{1} = 'fry';
exptAll{1} = 'u005_003'; % 

%%
slowXcorrFlag = set(G_handles.slowXcorr,'Value',0);  %lateral movement correction
fastXcorrFlag = set(G_handles.fastXcorr,'Value',1);  %lateral movement correction
OpticFlowFlag = set(G_handles.OpticFlowCorrection,'Value',1);  %lateral movement correction

frameRemovalFlag = set(G_handles.frameRemoval,'Value',1);  %lateral movement correction
deletionThresh_Z = set(G_handles.deletionThreshold,'string','0.5');

SVDremovalFlag = set(G_handles.PCremoval,'Value',0);  %lateral movement correction
N_PC = set(G_handles.nPC,'string','1');

DSflag = set(G_handles.DSflag,'Value',0);

%%

global f0m Tens


exdom = 1:length(exptAll);

for eid = 1:length(exdom)
    
    ex = exdom(eid);
    
    anim = animAll{ex};
    expt = exptAll{ex};
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    
    setGUIlabels
    
    %%%%%load 'cellS' for the motion information%%%%%%%
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
    load(tracepath,'cellS')
 
    
    set(G_handles.basesub,'Value',1);
    set(G_handles.bstart,'String','-300');  %Frame start in ms (to average)
    set(G_handles.bstop,'String','100'); %Frame stop in ms (to average)
    set(G_handles.epistart,'String','100'); %in msec as well
    set(G_handles.epistop,'String',num2str(getparam('stim_time')*1000+500)); %in msec as well
    
    processButton(1)
    
    
    
    %%%%%save 'f0m' and 'Tens'%%%%%%%
    traceroot = 'f:\2pF0preprocessed\';
    tracepath = [traceroot anim '_' expt(1:8) '_F0'];
    save(tracepath,'f0m','Tens')
    
end
