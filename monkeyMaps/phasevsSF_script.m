pF0
pRev

global G_RChandles G_handles PW NB TC

set(G_RChandles.kernelLength,'string','[-200 1000]');
set(G_RChandles.LPflag,'value',0);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);
set(G_RChandles.blankNorm,'value',0);



dataRoot = 'd:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';


%% ax3; u009_009

%Initialize structures with the first animal

clear dori dsf doriEuc dorisf ax dist doripair dsfpair dtcoripair dtcsfpair animID

global maskS cellS DM idExamp

%idExamp = [6 20 36];
idExamp = []
anim = 'ax3';
expt = 'u009_009_rigid';

%load expt
set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
Gsetdirectories

maskroot = 'C:\CellMasks\';
maskpath = [maskroot anim '_' expt(1:8)];
traceroot = 'C:\CellTraces\';
tracepath = [traceroot anim '_' expt(1:8) '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[21:65]')

Gkernelplots4;

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
animID = ones(1,length(TC.opref{1}));  %for bootstrap later

%close all
%% Get the rest of them

global maskS cellS

clear animS exptS
%animS{1} = 'ax3'; 
animS{2} = 'ax3'; 
animS{3} = 'ax2'; 
animS{4} = 'ax3';
animS{5} = 'ax3'; 
animS{6} = 'ax3'; 

%exptS{1} = 'u009_009'; %site3 (v1)
exptS{2} = 'u009_017'; %site3 (v1) Nice maps
exptS{3} = 'u000_108';  %Keith (v1)

exptS{4} = 'u009_027';  %site2 (v2?) Site 2; medial anterior. Some decent kernels
exptS{5} = 'u002_000'; %site2 (v2?) Some decent kernels
exptS{6} = 'u010_024'; %site4 (v2?) %Beautiful labeling, but only a few responsive cells

%ax3_10_21
%ax3_10_22

exdom = 2:2


idExampAll = [6 20 36; 12 28 57; 6 22 26; 23 31 50; 1 1 1; 1 1 1];

for exid = 1:length(exdom)
    
    %idExamp = idExampAll(exdom(exid),:);
    ex = exdom(exid);

    anim = animS{ex};
    expt = exptS{ex};

    %load expt
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories

    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
    
    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    set(G_RChandles.dropTrials,'string','[]')

    Gkernelplots4

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
    
    animID = [animID ex*ones(1,length(TC.opref{1}))];  %for bootstrap later
    
end

%%

%% Fit Gaussian
[sfpref2 id] = sort(sfpref);
F1F02 = F1F0(id);
id = find(~isnan(F1F02.*sfpref2));
F1F02 = F1F02(id);
sfpref2 = sfpref2(id);

sigdom = linspace(.1,3,100);
for i = 1:length(sigdom)
    ffit = 2*exp(-(sfpref2.^2)/(2*sigdom(i)^2));
    E(i) = norm(ffit-F1F02);
end
[dum id] = min(E);
sig = sigdom(id);
ffit = 2*exp(-(sfpref2.^2)/(2*sig^2));
figure,
subplot(1,2,1)
plot(sfpref2,ffit), hold on, plot(sfpref2,F1F02,'.k')
xlabel('SF preference'), ylabel('F1/F0')
title(['sig = ' num2str(sig) 'cyc/deg'])

%     [sfpref2 id] = sort(TC.sfpref2{c});
%     F1F02 = TC.F1F02{c}(id);
%     id = find(~isnan(F1F02.*sfpref2));
%     F1F02 = F1F02(id);
%     sfpref2 = sfpref2(id);
%     sfpref2 = [fliplr(-sfpref2) sfpref2];
%     F1F02 = [fliplr(F1F02) F1F02];
%     [param ffit varacc sigma] = Gaussfit2(sfpref2,F1F02);
%     figure,scatter(sfpref2,ffit), hold on, plot(sfpref2,F1F02,'.')


subplot(1,2,2), hist(F1F02,20), xlabel('F1/F0')
