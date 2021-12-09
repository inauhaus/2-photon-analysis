pF0

%initialize the Gui

global G_handles Analyzer cellS maskS ACQinfo twophDATADIR

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

%%
for eid = 1:length(exdom)
    
    ex = exdom(eid);
    
    anim = animAll{ex};
    expt = exptAll{ex};
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    setGUIlabels
   
   %%%%%%%%%%
    
    Tdelim = str2num(get(G_handles.maskTemplateDelim,'string'));
    fp = ACQinfo.linesPerFrame*ACQinfo.msPerLine/1000;
    predelayFrames = round(getparam('predelay')/fp);
    
    im = 0;
    
    fr = ACQinfo.SBInfo.resfreq/ACQinfo.SBInfo.sz(1); %frames/sec
    for i = 1:length(Tdelim)-1
        Twin = [Tdelim(i) Tdelim(i+1)];
        Fwin = round(Twin*fr)+1;
        dumdum = double(squeeze(sbxread(twophDATADIR,Fwin(1),1+Fwin(2)-Fwin(1))));
        dum{1} = dumdum(:,ACQinfo.unblanked,:);
        
        dum{1} = cleanTensor4(dum{1},Fwin,2); %This will downsample it in space and time.
        
        im = im + localXCorr3(dum{1},2)/(length(Tdelim)-1);
        %im = im + localXCorr4(dum{1},3);
    end
    
    %im = LocalZ(im,100,1);
    maskS.im{1} = im;
    
    dim = size(maskS.im{1});
    dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
    maskS.im{1} = interp1(1:dim(1),maskS.im{1},linspace(1,dim(1),dimI(1)));
    maskS.im{1} = interp1(1:dim(2),maskS.im{1}',linspace(1,dim(2),dimI(2)))';
    
    maskS.bw = zeros(size(maskS.im{1}));
    
    
    [xmicperpix ymicperpix] = getImResolution;
    xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
    ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
    figure(40),
    maskS.im{1} = (maskS.im{1} - mean(maskS.im{1}(:)))/std(maskS.im{1}(:));
    maskS.im{1}(find(maskS.im{1}<-2)) = -2;
    maskS.im{1}(find(maskS.im{1}>6)) = 6;
    imagesc(xdom,ydom,maskS.im{1}), colormap gray
    axis image
    hold off

    
    maskroot = 'C:\CellMasks\';
    maskpath = [traceroot anim '_' expt(1:8)];
    save(maskpath,'maskS')
    

end
