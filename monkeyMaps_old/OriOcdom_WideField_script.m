global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp flipeyebit SelectivityThresh

set(G_handles.datadir,'string','f:\2p_data\')
set(G_handles.analyzedir,'string','f:\2p_data\AnalyzerFiles\')

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

SelectivityThresh = .4;
%% ab8; u001_017 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run this if you want it to have the same ROI as the corresponding SF/Ocdom
%expt. i.e. get dprimemask and use it as an input 
idExamp = []; 
anim = 'ab8'; dpthresh = 1;   flipeyebit = 0;
[dum1 dum2 maskdum dum3] = CombineOcdomSFexpt(anim, dpthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idExamp = [90 190; 200 217; 200 140; 216 98];  %left hemi

anim = 'ab8';
expt = 'u001_017';

eno = 1;

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1; flipeyebit = 0;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] = OriOcdom_WideField(dthresh,maskdum);

set(G_handles.Lwidth,'string','.eno');
hh = makeMapFilter;
plotMapExamples_oriocdom(hh) %shows ori and od curves


if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%% ab9; u000_102 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run this if you want it to have the same ROI as the corresponding SF/Ocdom
%expt. i.e. get dprimemask and use it as an input 
idExamp = []; 
anim = 'ab9'; dpthresh = 3;   flipeyebit = 1;
[dum1 dum2 maskdum dum3] = CombineOcdomSFexpt(anim, dpthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eno = 2;

anim = 'ab9';
expt = 'u000_102';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3; flipeyebit = 1;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] = OriOcdom_WideField(dthresh,maskdum);

if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%% ac0; u000_056 This one is kinda shitty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run this if you want it to have the same ROI as the corresponding SF/Ocdom
%expt. i.e. get dprimemask and use it as an input 
% idExamp = []; 
% anim = 'ac0'; dpthresh = 1;   flipeyebit = 1;
% [dum1 dum2 maskdum dum3] = CombineOcdomSFexpt(anim, dpthresh);

anim = 'ac0';  %This one is really noisy, but still consistent (may want to  smooth more)
expt = 'u000_057';  %ORI4,2BW,5SF
%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels
load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
dthresh = 1; flipeyebit = 1;  %left hemisphere;
[dum1 dum2 maskdum dum3] = SfOcdom_WideField(dthresh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eno = 3;

anim = 'ac0';
expt = 'u000_056';  %ori8,2BW
%N.B. I can't use the orientation map from 57 because there are only 4
%orientations

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1; flipeyebit = 1;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] = OriOcdom_WideField(dthresh,maskdum);


if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%% The rest are from that also show different SFs 

%%
% anim = 'ac0';  %This one is garbage
% expt = 'u000_063';  
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% setGUIlabels
% 
% try
%     load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])    
% catch
%     load(['e:\Beta\f0images\f0_' anim '_' expt(2:end)])
% end
% 
% dthresh = .1
% [dum ocdommap_raw dprimemask dum2 ocdommap] = OriOcdom_WideField(dthresh);



%%

anim = 'ac1';
expt = 'u000_115';

eno = 4;

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])    
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1; flipeyebit = 1;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] =  OriOcdom_WideField(dthresh);

if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%%

anim = 'ac1';
expt = 'u001_011';

eno = 5;

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1; flipeyebit = 0;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] = OriOcdom_WideField(dthresh);

if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%%
eno = 6;

anim = 'ac1';
expt = 'u002_019';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1;  flipeyebit = 0;
[dum ocdommap_raw{eno} dprimemask{eno} dum2 ocdommap{eno} ORIODdangROI{eno} ORIODdangROISEL{eno}] = OriOcdom_WideField(dthresh);

if CCbit
    %[oriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dprimemask{eno},1);
    [odCC{eno} odCCdom{eno}] = circcorr2D(ocdommap{eno},dprimemask{eno},0);
end

%% Combine experiments


%% Intersection histograms
intangleROIAll = []; intangleROISELAll = []; TotalA = 0;
for i = 1:length(ORIODdangROI)    
   intangleROIAll = [intangleROIAll; ORIODdangROI{i}];    
   intangleROISELAll = [intangleROISELAll; ORIODdangROISEL{i}];  
   
   TotalA = length(find(dprimemask{i})) + TotalA;
end

%length(intangleROISELAll)/TotalA   %Percentage of ROI area that had "reliable gradients" in both OD and SF maps
plotGradIntersection(intangleROIAll,intangleROISELAll)

%%

if CCbit
    minrange = find(oriCCdom{1}<600);
    figure
   for i = 1:length(sfCC)
      
       subplot(length(sfCC),2,(i-1)*2+1)
       plot(oriCCdom{i},oriCC{i})
       xlim([0 650]), ylim([-.6 1])
       title('Ori acorr'), xlabel('um')
       
       subplot(length(odCC),2,i*2)
       plot(odCCdom{i},odCC{i})
       xlim([0 650]), ylim([-.6 1])
       title('Ocdom acorr'), xlabel('um')
       
       [dum orimin(i)] = min(oriCC{i}(minrange));
       orimin(i) = oriCCdom{1}(orimin(i));
       
       [dum odmin(i)] = min(odCC{i}(minrange));
       odmin(i) = odCCdom{1}(odmin(i));
       
   end    
   
   figure,scatter(odmin,orimin,'k')
   hold on, plot([0 600],[0 600],'r')
   xlabel('ocular dominance half period'), ylabel('orientation half period')
    
end