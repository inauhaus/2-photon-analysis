global G_RChandles G_handles PW NB TC

%% ab2; u000_017

clear dori dpos doriNorm dposNorm doripos ax dist sizesum

global maskS cellS

anim = 'ab2';
expt = 'u000_017';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3

%The first experiment analyzed initializes these structures
%Get pairwise info
for i = 1:length(PW.dori)  %loop through each distance
    dori{1}{i} = PW.dori{i}(:);
    dpos{1}{i} = PW.dpos{i}(:);
    sizesum{1}{i} = PW.sizesum{i}(:);
    doriNorm{1}{i} = PW.doriNorm{i}(:);
    dposNorm{1}{i} = PW.dposNorm{i}(:);
    ax{1}{i} = PW.ax{i}(:);
    dist{1}{i} = PW.dist{i}(:);
end


close all
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
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3

for i = 1:length(PW.dori)  %loop through each distance
    dori{2}{i} = PW.dori{i}(:);
    dpos{2}{i} = PW.dpos{i}(:);
    sizesum{2}{i} = PW.sizesum{i}(:);
    doriNorm{2}{i} = PW.doriNorm{i}(:);
    dposNorm{2}{i} = PW.dposNorm{i}(:);
    ax{2}{i} = PW.ax{i}(:);
    dist{2}{i} = PW.dist{i}(:);
end


close all
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
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3

for i = 1:length(PW.dori)  %loop through each distance
    dori{3}{i} = PW.dori{i}(:);
    dpos{3}{i} = PW.dpos{i}(:);
    sizesum{3}{i} = PW.sizesum{i}(:);
    doriNorm{3}{i} = PW.doriNorm{i}(:);
    dposNorm{3}{i} = PW.dposNorm{i}(:);
    ax{3}{i} = PW.ax{i}(:);
    dist{3}{i} = PW.dist{i}(:);
end


close all



%% aa9; u001_031

global maskS cellS

anim = 'aa9';
expt = 'u001_031';
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

%load analyzer file
anim = 'aa9';
expt = 'u001_031';
dir = get(G_handles.analyzedir,'String'); %partial path for analyzer file
setAnalyzerDirectory([dir anim '\' ])
loadAnalyzer(expt)

%get the ACQinfo
expt = 'u001_027';
dir = get(G_handles.datadir,'String'); %partial path for .tiff files 
twophDATADIR = [dir anim '\' expt '\'];  %Make path for .tiff files
AUE = [anim '_' expt]; %loaded unit experiment 'u000_000'
setacqinfo(1)  %Set the global ACQinfo structure (contains trial info

Grandposplots3

for i = 1:length(PW.dori)  %loop through each distance
    dori{4}{i} = PW.dori{i}(:);
    dpos{4}{i} = PW.dpos{i}(:);
    sizesum{4}{i} = PW.sizesum{i}(:);
    doriNorm{4}{i} = PW.doriNorm{i}(:);
    dposNorm{4}{i} = PW.dposNorm{i}(:);
    ax{4}{i} = PW.ax{i}(:);
    dist{4}{i} = PW.dist{i}(:);
end

close all

%% ab3; u002_017 (near edge)

global maskS cellS

anim = 'ab3';
expt = 'u002_017';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3

for i = 1:length(PW.dori)  %loop through each distance
    dori{5}{i} = PW.dori{i}(:);
    dpos{5}{i} = PW.dpos{i}(:);
    sizesum{5}{i} = PW.sizesum{i}(:);
    doriNorm{5}{i} = PW.doriNorm{i}(:);
    dposNorm{5}{i} = PW.dposNorm{i}(:);
    ax{5}{i} = PW.ax{i}(:);
    dist{5}{i} = PW.dist{i}(:);
end

close all






% %% ab4; u000_027
% 
% global maskS cellS
% 
% anim = 'ab4';
% expt = 'u000_027';
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% 
% maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
% maskpath = [maskroot anim '_' expt];
% traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
% tracepath = [traceroot anim '_' expt '_cellS'];
% 
% load(maskpath,'maskS') 
% load(tracepath,'cellS') 
% 
% set(G_RChandles.dropTrials,'string','[]')
% 
% Grandposplots3
% 
% for i = 1:length(PW.dori)  %loop through each distance
%     dori{6}{i} = PW.dori{i}(:);
%     dpos{6}{i} = PW.dpos{i}(:);
%     sizesum{6}{i} = PW.sizesum{i}(:);
%     doriNorm{6}{i} = PW.doriNorm{i}(:);
%     dposNorm{6}{i} = PW.dposNorm{i}(:);
%     ax{6}{i} = PW.ax{i}(:);
%     dist{6}{i} = PW.dist{i}(:);
% end
% 
% close all

%% ab3; u002_028 (on edge of stimulus)
% 
% global maskS cellS
% 
% anim = 'ab3';
% expt = 'u002_028';
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% 
% maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
% maskpath = [maskroot anim '_' expt];
% traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\';
% tracepath = [traceroot anim '_' expt '_cellS'];
% 
% load(maskpath,'maskS') 
% load(tracepath,'cellS') 
% 
% set(G_RChandles.dropTrials,'string','[]')
% 
% Grandposplots3
% 
% for i = 1:length(PW.dori)  %loop through each distance
%     dori{7}{i} = PW.dori{i}(:);
%     dpos{7}{i} = PW.dpos{i}(:);
%     sizesum{7}{i} = PW.sizesum{i}(:);
%     doriNorm{7}{i} = PW.doriNorm{i}(:);
%     dposNorm{7}{i} = PW.dposNorm{i}(:);
%     ax{7}{i} = PW.ax{i}(:);
%     dist{7}{i} = PW.dist{i}(:);
% end
% 
% close all
% 
