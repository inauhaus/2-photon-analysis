%This OriSf script (version 2) was made to compile the ORI/SF experiments that also have
%Ocular dominance (i.e. for the OD vs. SF paper)
%%

global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp flipeyebit SelectivityThresh

set(G_handles.datadir,'string','f:\2p_data\')
set(G_handles.analyzedir,'string','f:\2p_data\AnalyzerFiles\')

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

SelectivityThresh = .4;

%%%%%%%First, left hemisphere%%%%%%%%%%%

idExamp = [100 217; 211 92; 156 140; 138 191];  %left hemi


CCbit = 0;
%% ab8; u001_007   (~230um deep; middle)
eno = 1;

anim = 'ab8';
expt = 'u001_007';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

idExamp = [90 210; 200 217; 216 98; 200 140];  %left hemi
dthresh = 1;
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.1');
hh = makeMapFilter;
plotMapExamples(hh)


if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end

close all
%% ab9; u000_084 (bigger spot, thats decent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run this if you want it to have the same ROI as the corresponding SF/Ocdom
%expt. i.e. get dprimemask and use it as an input
idExamp = [];
anim = 'ab9'; dpthresh = 3;   flipeyebit = 1;
[dum1 dum2 maskdum dum3] = CombineOcdomSFexpt(anim, dpthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eno = 2;

anim = 'ab9';
expt = 'u000_084';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)]) %N.B.  A trial from experiment 80 (aborted) was averaged in with f0m and f0m_var.

dthresh = 3;
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh,maskdum);

if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end

%close all
%% ac0; u000_047  (ok)

eno = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run this if you want it to have the same ROI as the corresponding SF/Ocdom
%expt. i.e. get dprimemask and use it as an input
anim = 'ac0';  %This one is really noisy, but still consistent (may want to  smooth more)
expt = 'u000_057';  %ORI4,2BW,5SF
%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels
load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
dthresh = 1; flipeyebit = 1;
[dum1 dum2 maskdum dum3] = SfOcdom_WideField(dthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anim = 'ac0';
expt = 'u000_047';   %ORI8/SF5
%N.B. Can't use expt 57 here because it only has 4 oris

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1;
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh,maskdum);

if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end
%close all

%%

eno = 4;

anim = 'ac1';
expt = 'u000_115';  %ORI8/SF5/OD2

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

dthresh = 1;
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh);

if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end
close all

%%

eno = 5;

anim = 'ac1';
expt = 'u001_011';  %ORI8/SF5/OD2

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

dthresh = 1;

idExamp = [];
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh);

if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end

%%

eno = 6;

anim = 'ac1';
expt = 'u002_019';  %ORI8/SF5/OD2

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

dthresh = 1;
[oriang{eno} orimag{eno} sfpref{eno} sfmag{eno} dpmask{eno} ORISFdangROI{eno} ORISFdangROISEL{eno}] = OriSf_WideField(dthresh);

if CCbit
    [oriCC{eno} imoriCC{eno} oriCCdom{eno}] = circcorr2D(oriang{eno},dpmask{eno},1);
    [sfCC{eno} imsfCC{eno} sfCCdom{2}] = circcorr2D(log2(sfpref{eno}),dpmask{eno},0);
end

close all

%% Combine experiments

%% Intersection histograms
intangleROIAll = []; intangleROISELAll = []; TotalA = 0;
for i = 1:length(ORISFdangROI)    
   intangleROIAll = [intangleROIAll; ORISFdangROI{i}];    
   intangleROISELAll = [intangleROISELAll; ORISFdangROISEL{i}];  
   
   TotalA = length(find(dprimemask{i})) + TotalA;
end

plotGradIntersection(intangleROIAll,intangleROISELAll)
%length(intangleROISELAll)/TotalA   %Percentage of ROI area that had "reliable gradients" in both OD and SF maps

%%

if CCbit
    figure
    for i = 1:length(sfCC)

        subplot(length(sfCC),2,(i-1)*2+1)
        plot(oriCCdom{i},oriCC{i})
        xlim([0 650]), ylim([-.3 1])
        title('Ori acorr'), xlabel('um')

        subplot(length(sfCC),2,i*2)
        plot(sfCCdom{i},sfCC{i})
        xlim([0 650]), ylim([-.3 1])
        title('Sf acorr'), xlabel('um')

    end

end