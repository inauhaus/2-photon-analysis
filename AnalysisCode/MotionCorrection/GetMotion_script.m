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
exptAll{id} = 'u005_003';

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


%%%%%%%%%%%%%%%
id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u011_002'; %Spiderman

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u011_003'; %texture v noise

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u011_004'; %random gratings

id = id+1;
animAll{id} = 'fry';
exptAll{id} = 'u011_005'; %periodic retinotopy


%%
exdom = 1:length(exptAll);
exdom = 14:26
%%
for eid = 1:length(exdom)
    
    ex = exdom(eid);
    
    anim = animAll{ex};
    expt = exptAll{ex};
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    
    
    %%%%%%%%%%

    
    T = 3; %seconds
    nT = 30; %N windows
    
    [bestWindow maskS.anatomyTemplate] = findTempWindow(T,nT);
    
    mi = prctile(maskS.anatomyTemplate(:),1);
    ma = prctile(maskS.anatomyTemplate(:),99.8);
    
     figure,imagesc(maskS.anatomyTemplate(:,:),[mi ma]), colormap gray
     title('Anatomy Template')
     axis image
     drawnow

    cellS.anatomyTemplate = maskS.anatomyTemplate;  %Better to save it in the location that has the motion correction, so I know they match up.
      
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    Win = 10;  %seconds
    experimentMotionAnalyzer2(Win)
    
    %%%%%save 'cellS' w/o the traces, only motion information%%%%%%%
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
    save(tracepath,'cellS')
 
    
end
