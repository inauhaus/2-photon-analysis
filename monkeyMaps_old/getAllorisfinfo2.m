global G_RChandles G_handles PW NB TC

%% ab2; u000_014  

%Initialize structures with the first animal

clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID

global maskS cellS

anim = 'ab2';
expt = 'u000_014';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

[maporidiff(1) mapslopemag(1)] = Gkernelplots4

%Initialize pairwise info
fnames = fieldnames(PW);
for i = 1:length(fnames)
    dum = eval(['PW.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    else
        eval([fnames{i} '= dum;']);
    end
end
animID_d = ones(1,length(PW.dori{1}));  %for bootstrap later (identifies each pair with an animal)

%Initialize neighborhood info
fnames = fieldnames(NB);
for i = 1:length(fnames)
    dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    end
end

%Initialize single-cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    else
        eval([fnames{i} '= dum;']);
    end
end
animID = ones(1,length(TC.opref{1}));  %for bootstrap later

close all
%% Get the rest of them

global maskS cellS

clear animS exptS
%animS{1} = 'ab2';  %160
animS{2} = 'ab2'; 
animS{3} = 'ab2'; 
animS{4} = 'ab2';
animS{5} = 'ab3'; 
animS{6} = 'ab3'; 
animS{7} = 'aa9'; 
animS{8} = 'ab4';

animS{9} = 'ab3'; 
animS{10} = 'ab3';

%exptS{1} = 'u000_014'; %Small spot, but nice orthogonal linear zones;  160 um deep 
exptS{2} = 'u000_068'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u000_087'; %small region. linear zone and orth sfreq
exptS{4} = 'u001_004'; %fast linear zone, orthogonal sfreq map; 160 um deep
exptS{5} = 'u002_014'; %isodomain, nice spat freq map ; 100 um deep (LMS is 2_13; randpos is 2_17)
exptS{6} = 'u002_041'; %isodomain, awesome spat freq map (above/below 2_14); 140 um deep
exptS{7} = 'u001_027';  %pinwheel, low spat freq; 175 um deep
exptS{8} = 'u001_003'; %ori fracture, low spat freq; 170 um deep

exptS{9} = 'u002_027'; %isodomain, sp freq linear zone ; 200 um deep
exptS{10} = 'u004_016';

for ex = 2:10

    anim = animS{ex};
    expt = exptS{ex};

    %load expt
    set(G_handles.loadana,'string',anim)
    set(G_handles.loadexp,'string',expt)
    Gsetdirectories

    maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
    maskpath = [maskroot anim '_' expt];
    traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
    tracepath = [traceroot anim '_' expt '_cellS'];

    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    if strcmp(anim,'ab2') & strcmp(expt,'u000_087')
        set(G_RChandles.dropTrials,'string','[14:30]')
    else
        set(G_RChandles.dropTrials,'string','[]')
    end

    [maporidiff(ex) mapslopemag(ex)] = Gkernelplots4
    close all
    
    %Accumulate pairwise info
    fnames = fieldnames(PW);
    for i = 1:length(fnames)
        dum = eval(['PW.' fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:)];']);
        else
            eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
    animID_d = [animID_d ex*ones(1,length(PW.dori{1}))];  %for bootstrap later  (identifies each pair with an animal)

    %Accumulate neighborhood info
    fnames = fieldnames(NB);
    for i = 1:length(fnames)
        dum = eval(['NB.' fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:);]']);
        end
    end

    %Accumulate single cell info
    fnames = fieldnames(TC);
    for i = 1:length(fnames)
        dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:)];']);
        else
            eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
    
    animID = [animID ex*ones(1,length(TC.opref{1}))];  %for bootstrap later
    
end



%%

% 
% 
clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID idExamp

global maskS cellS idExamp

anim = 'ab2';
expt = 'u000_014';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

%idExamp = [6 27]
%idExamp = [79 59 33]
%idExamp = [6 34 54 79];  %ab3 002_041, ab2  000_068
%idExamp = [9 54 78 83];  %aa9 001_027
idExamp = [3 31 54];  %ab2 000_014

Gkernelplots4
% 


% 

% global maskS cellS idExamp
% 
% anim = 'ab2';
% expt = 'u000_014';
% 
% %load expt
% set(G_handles.loadana,'string',anim)
% set(G_handles.loadexp,'string',expt)
% Gsetdirectories
% 
% maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
% maskpath = [maskroot anim '_' expt(1:8)];
% traceroot = 'C:\2ph_code\Beta\cellTraces\monkey\new\';
% tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
% 
% load(maskpath,'maskS') 
% load(tracepath,'cellS') 
% 
% set(G_RChandles.dropTrials,'string','[]')
% 
% idExamp = [];
% 
% Gkernelplots4
% 
% 
% 
