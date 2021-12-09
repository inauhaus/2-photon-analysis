function getAlloriposinfo

%% ax3; u009_015

%Initialize structures with the first animal

clear dori dpos doriNorm dposNorm doripos ax dist sizesum animID

global maskS cellS idExamp

idExamp = [5 8 10 11 12 20 21 25];

anim = 'ax3';
expt = 'u009_015';
c = 1;
kElem = 1;

%load expt
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) '_rigid'])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

Grandposplots3;

%Initialize pairwise info
fnames = fieldnames(PW);
for i = 1:length(fnames)
    dum = eval(['PW.' fnames{i} '{c}{kElem};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
    else
        eval([fnames{i} '= dum;']);
    end
end
animID_d = ones(1,length(PW.dori{c}{kElem}));  %for bootstrap later (identifies each pair with an animal)

%Initialize single-cell info
fnames = fieldnames(TC);
for i = 1:length(fnames)
    dum = eval(['TC.' fnames{i} '{c}{kElem};']); dim = size(dum); 
    if min(dim) == 1
        eval([fnames{i} '= dum(:);']);
%     else
%         eval([fnames{i} '= dum;']);
    end
end
animID = ones(1,length(TC.xpos{c}{kElem}));  %for bootstrap later

%PI{i} = PIs;
 
%close all
%% Get the rest of them

global maskS cellS

clear animS exptS
%animS{1} = 'ax3'; %Zeus 
animS{2} = 'ax3'; 
animS{3} = 'ax3'; 
animS{4} = 'ax3';

animS{5} = 'ax2'; %Keith
animS{6} = 'ax2';
animS{7} = 'ax2'; %Keith
animS{8} = 'ax2';

%animS{9} = 'ax3';



%exptS{1} = 'u009_015'; V1 (Zeus)
exptS{2} = 'u009_010'; %V1 (Zeus) 19 good kernels
exptS{3} = 'u009_021';  %V2 (Zeus) 
exptS{4} = 'u009_022';  %V2 (Zeus) about 6 good kernels

exptS{5} = 'u000_102'; %Keith V1 (these stink)
exptS{6} = 'u000_106'; %Keith V1 about 5 kernels.  Very few cells in mask
exptS{7} = 'u000_092'; %Keith. V1 Ok.  Off center though.  On off clustering ok
exptS{8} = 'u000_095'; %Keith only one ori.  Very clear On off clustering

%exptS{9} = 'u001_035';  %V2

clear idExampAll
idExampAll{1} = [5 10 11 25];
idExampAll{2} = [10 11 12 14];
idExampAll{3} = [1 16 26 29 31];
idExampAll{4} = [15 17 18];
idExampAll{6} = [1 3 4 5 7];
idExampAll{7} = [13 15 22 28 29 42]
idExampAll{8} = []


exdom = [3 4 7]; 
%exdom = [2 5 6] 
exdom = 6


for exid = 1:length(exdom)
    try
        idExamp = idExampAll{exdom(exid)};
    catch
        idExamp = [];
    end
    
    
    ex = exdom(exid);

    anim = animS{ex};
    expt = exptS{ex};
    
    %load expt
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    
    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
    
    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    set(G_RChandles.dropTrials,'string','[]')
    if strcmp(expt,'u000_092') & strcmp(anim,'ax2')
       set(G_RChandles.dropTrials,'string','[24:930]');
    end
    
    
    Grandposplots3
    
    %Accumulate pairwise info
    fnames = fieldnames(PW);
    for i = 1:length(fnames)
        dum = eval(['PW.' fnames{i} '{c}{kElem};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:)];']);
        else
            eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
    animID_d = [animID_d ex*ones(1,length(PW.dori{c}{kElem}))];  %for bootstrap later  (identifies each pair with an animal)

    %Accumulate single cell info
    fnames = fieldnames(TC);
    for i = 1:length(fnames)
        dum = eval(['TC.' fnames{i} '{c}{kElem};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '= [' fnames{i} '; dum(:)];']);
%         else
%             eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
    
    animID = [animID ex*ones(1,length(TC.xpos{c}{kElem}))];  %for bootstrap later
    
    %PI{ex} = PIs;
    
end

%% 
