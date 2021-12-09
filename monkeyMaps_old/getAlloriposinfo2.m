global G_RChandles G_handles PW TC PI PIs

%% ab2; u000_017

%Initialize structures with the first animal

clear dori dpos doriNorm dposNorm doripos ax dist sizesum animID

global maskS cellS

anim = 'ab2';
expt = 'u000_017';

% anim = 'ab2';
% expt = 'u001_006';

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

Grandposplots3;

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

%Initialize single-cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
%     else
%         eval([fnames{i} '= dum;']);
    end
end
animID = ones(1,length(TC.xpos{1}));  %for bootstrap later

PI{i} = PIs;

%close all
%% Get the rest of them

global maskS cellS

clear animS exptS
%animS{1} = 'ab2'; 
animS{2} = 'ab2'; 
animS{3} = 'ab2'; 
animS{4} = 'aa9';
animS{5} = 'ab3'; 
animS{6} = 'ab8';
animS{7} = 'ab4'; 

animS{8} = 'ab4'; 
animS{9} = 'ab3'; 


%exptS{1} = 'u000_017'; %Small spot, but nice orthogonal linear zones
exptS{2} = 'u000_070'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u001_006';  %Pinwheel (used to look like a linear zone with ori/sf... maybe shifted ROI)
exptS{4} = 'u001_031';  %pinwheel, low spat freq
exptS{5} = 'u002_017'; %isodomain, nice spat freq map
exptS{6} = 'u001_035';  %linear zone with not much transition, nice On/Off map
exptS{7} = 'u001_050'; %pinwheel


exptS{8} = 'u000_027'; 
exptS{9} = 'u002_028';  %poorly centered

exdom = [3 4 7]; %pinwheels
%exdom = [2 5 6] %non pinwheels
exdom = 4:4

for exid = 1:length(exdom)
    ex = exdom(exid);

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

    Grandposplots3

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

    %Accumulate single cell info
    fnames = fieldnames(TC);
    for i = 1:length(fnames)
        dum = eval(['TC.' fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:)];']);
%         else
%             eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
    
    animID = [animID ex*ones(1,length(TC.xpos{1}))];  %for bootstrap later
    
    PI{ex} = PIs;
    
end

%% 

%Initialize structures with the first animal

clear dori dpos doriNorm dposNorm doripos ax dist sizesum animID

idExamp = 55;  %ab2 0_17; this is the same as 54 for 0_14 mask


global maskS cellS

anim = 'zc3';
expt = 'u000_011';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories

maskroot = 'C:\2ph_code\Beta\Masks\mouse\';
maskpath = [maskroot anim '_' expt];
traceroot = 'C:\2ph_code\Beta\cellTraces\mouse\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3

