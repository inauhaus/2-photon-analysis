global G_RChandles G_handles PW NB TC

%% ab2; u000_011

%Initialize structures with the first animal

clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID

global maskS cellS DM

anim = 'ab2';
expt = 'u000_011';

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

if length(DM.colordom) == 3
    colordomdum = 1:4;
end

%Initialize pairwise info
fnames = fieldnames(PW);
for c = 1:length(colordomdum)
    for i = 1:length(fnames)
        dum = eval(['PW.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i}  '{' num2str(c) '} = dum(:);']);
        else
            eval([fnames{i}  '{' num2str(c) '} = dum;']);
        end
    end
end
animID_d = ones(1,length(PW.dori{1}));  %for bootstrap later (identifies each pair with an animal)

%Initialize neighborhood info
% fnames = fieldnames(NB);
% for c = 1:length(DM.colordom)
%     for i = 1:length(fnames)
%         dum = eval(['NB.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
%         if min(dim) == 1
%             eval([fnames{i} '{' num2str(c) '}= dum(:);']);
%         end
%     end
% end

%Initialize single-cell info
fnames = fieldnames(TC);
for c = 1:length(colordomdum)
    for i = 1:length(fnames)
        dum = eval(['TC.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
        if min(dim) == 1
            eval([fnames{i} '{' num2str(c) '} = dum(:);']);
        else
            eval([fnames{i} '{' num2str(c) '} = dum;']);
        end
    end
end
animID = ones(1,length(TC.opref{1}));  %for bootstrap later

close all
%% Get the rest of them

global maskS cellS

clear animS exptS
%animS{1} = 'ab2'; 
animS{2} = 'ab2'; 
animS{3} = 'ab2'; 
animS{4} = 'aa9';
animS{5} = 'ab3'; 
animS{6} = 'ab2'; 
animS{7} = 'ab3'; 



%exptS{1} = 'u000_011'; %Small spot, but nice orthogonal linear zones
exptS{2} = 'u000_057'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u000_089'; %small region. linear zone and orth sfreq
exptS{4} = 'u001_026';  %pinwheel, low spat freq
exptS{5} = 'u002_013'; %isodomain, nice spat freq map
exptS{6} = 'u002_018'; %fast linear zone, orthogonal sfreq map
exptS{7} = 'u002_042'; %isodomain, awesome spat freq map (above/below 2_14)


for ex = 2:6

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
    
    set(G_RChandles.dropTrials,'string','[]')

    [maporidiff(ex) mapslopemag(ex)] = Gkernelplots4
    close all
    
    %Accumulate pairwise info
    fnames = fieldnames(PW);
    for c = 1:length(colordomdum)
        for i = 1:length(fnames)
            dum = eval(['PW.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
            if min(dim) == 1
                eval([fnames{i} '{' num2str(c) '} = [' fnames{i} '{' num2str(c) '}; dum(:)];']);
            else
                eval([fnames{i} '{' num2str(c) '} = [' fnames{i} '{' num2str(c) '}; dum];']);
            end
        end
    end
    animID_d = [animID_d ex*ones(1,length(PW.dori{1}))];  %for bootstrap later  (identifies each pair with an animal)

    %Accumulate neighborhood info
%     fnames = fieldnames(NB);
%     for c = 1:length(colordomdum)
%         for i = 1:length(fnames)
%             dum = eval(['NB.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
%             if min(dim) == 1
%                 eval([fnames{i} '{' num2str(c) '} = [' fnames{i} '{' num2str(c) '}; dum(:)];']);
%             end
%         end
%     end

    %Accumulate single cell info
    fnames = fieldnames(TC);
    for c = 1:length(colordomdum)
        for i = 1:length(fnames)
            dum = eval(['TC.' fnames{i} '{' num2str(c) '};']); dim = size(dum);
            if min(dim) == 1
                eval([fnames{i} '{' num2str(c) '} = [' fnames{i} '{' num2str(c) '}; dum(:)];']);
            else
                eval([fnames{i} '{' num2str(c) '} = [' fnames{i} '{' num2str(c) '}; dum];']);
            end
        end
    end
    
    animID = [animID ex*ones(1,length(TC.opref{1}))];  %for bootstrap later
    
end



%%

% 
% 
clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID idExamp

global maskS cellS idExamp

anim = 'ab3';
expt = 'u002_013';

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

idExamp = [44];  %ab2 0_11  (same as 54 for 0_14)

getTCfromRevCorr3

%Gkernelplots4
% 

