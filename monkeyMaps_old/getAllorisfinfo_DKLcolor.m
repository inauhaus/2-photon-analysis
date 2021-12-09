global G_RChandles G_handles PW NB TC

%% ab2; u000_015

%Initialize structures with the first animal

clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID

global maskS cellS DM

anim = 'ab3';
expt = 'u002_015';

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

%close all
%% Get the rest of them

global maskS cellS

clear animS exptS

%animS{1} = 'ab3'; 
animS{2} = 'ab3'; 
animS{3} = 'ab3';
animS{4} = 'ab4'; 
animS{5} = 'ab4'; 
animS{6} = 'ab4'; 
animS{7} = 'ab2'; 
animS{8} = 'ab4'; 

% animS{8} = 'ab2';
% animS{9} = 'ab3';
% animS{10} = 'ab3'; 



%exptS{1} = 'u002_015'; %
exptS{2} = 'u002_039'; %
exptS{3} = 'u004_017'; %
exptS{4} = 'u000_025'; %this one actually has a some ok tuning curves
exptS{5} = 'u001_047'; %
exptS{6} = 'u001_051'; %Trial 25 out of focus
exptS{7} = 'u000_062'; %really unstable after trial 30 (so aborted)

exptS{8} = 'u001_032'; %really unstable after trial 30 (so aborted)

% exptS{8} = 'u000_062'; %15 more trials in 63 after replacing the water
% exptS{9} = 'u002_054';  %14 trials (movement in trial 8).  Maybe 2_55 is better
% exptS{10} = 'u005_014'; %only one orientation (maybe 5_15 is better?)


for ex = 2:7

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
    if strcmp(expt,'u000_062')
        set(G_RChandles.dropTrials,'string','[28:40]')
    end

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

