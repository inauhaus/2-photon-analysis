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
%traceroot = 'e:\Beta\cellTraces\monkey\new2\';
traceroot = 'e:\Beta\cellTraces\monkey\noZ\';
%traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

[maporidiff(1) mapslopemag(1)] = Gkernelplots4

%%%%%%%%%
xshift = 6; yshift = 7;
prcoverlap = getcontroloverlap(xshift,yshift);
id = find(prcoverlap > .01);
TC.kernE{1}(id) = NaN;
%%%%%%%%%

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
animS{1} = 'ab2'; 
animS{2} = 'ab2'; 
animS{3} = 'ab2';
animS{4} = 'ab3'; 
animS{5} = 'ab3'; 
animS{6} = 'aa9'; 
animS{7} = 'ab4';

animS{8} = 'ab2'; 
animS{9} = 'ab3'; 
animS{10} = 'ab3';

exptS{1} = 'u000_014'; %Small spot, but nice orthogonal linear zones
exptS{2} = 'u000_068'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u001_004'; %fast linear zone, orthogonal sfreq map
exptS{4} = 'u002_014'; %isodomain, nice spat freq map
exptS{5} = 'u002_041'; %isodomain, awesome spat freq map (above/below 2_14)
exptS{6} = 'u001_027';  %pinwheel, low spat freq
exptS{7} = 'u001_003'; %ori fracture, low spat freq

exptS{8} = 'u000_087'; %small region. linear zone and orth sfreq
exptS{9} = 'u002_027'; %isodomain, sp freq linear zone
exptS{10} = 'u004_016';

for ex = 2:2

    anim = animS{ex};
    expt = exptS{ex};

    %load expt
    set(G_handles.loadana,'string',anim)
    set(G_handles.loadexp,'string',expt)
    Gsetdirectories

    maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
    maskpath = [maskroot anim '_' expt];
     %traceroot = 'e:\Beta\cellTraces\monkey\new2\';
    traceroot = 'e:\Beta\cellTraces\monkey\noZ\';
    %traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
    tracepath = [traceroot anim '_' expt '_cellS'];

    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    if strcmp(anim,'ab2') & strcmp(expt,'u000_087')
        set(G_RChandles.dropTrials,'string','[14:30]')
    else
        set(G_RChandles.dropTrials,'string','[]')
    end

    [maporidiff(ex) mapslopemag(ex)] = Gkernelplots4
    
    %%%%%%%%%
    xshift = 6; yshift = 7;
    prcoverlap = getcontroloverlap(xshift,yshift);
    id = find(prcoverlap > .01);
    TC.kernE{1}(id) = NaN;
    %%%%%%%%%
    
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

SNRactual = kernE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now the control (shifted mask)

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
%traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
traceroot = 'i:\Beta\cellTraces\monkey\noZ_control\';
tracepath = [traceroot anim '_' expt '_cellS'];

load(maskpath,'maskS') 
load(tracepath,'cellS') 

set(G_RChandles.dropTrials,'string','[]')

[maporidiff(1) mapslopemag(1)] = Gkernelplots4

%%%%%%%%%
xshift = 6; yshift = 7;
prcoverlap = getcontroloverlap(xshift,yshift);
id = find(prcoverlap > .01);
TC.kernE{1}(id) = NaN;
%%%%%%%%%

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
animS{1} = 'ab2'; 
animS{2} = 'ab2'; 
animS{3} = 'ab2';
animS{4} = 'ab3'; 
animS{5} = 'ab3'; 
animS{6} = 'aa9'; 
animS{7} = 'ab4';

animS{8} = 'ab2'; 
animS{9} = 'ab3'; 
animS{10} = 'ab3';

exptS{1} = 'u000_014'; %Small spot, but nice orthogonal linear zones
exptS{2} = 'u000_068'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u001_004'; %fast linear zone, orthogonal sfreq map
exptS{4} = 'u002_014'; %isodomain, nice spat freq map
exptS{5} = 'u002_041'; %isodomain, awesome spat freq map (above/below 2_14)
exptS{6} = 'u001_027';  %pinwheel, low spat freq
exptS{7} = 'u001_003'; %ori fracture, low spat freq

exptS{8} = 'u000_087'; %small region. linear zone and orth sfreq
exptS{9} = 'u002_027'; %isodomain, sp freq linear zone
exptS{10} = 'u004_016';

for ex = 2:7

    anim = animS{ex};
    expt = exptS{ex};

    %load expt
    set(G_handles.loadana,'string',anim)
    set(G_handles.loadexp,'string',expt)
    Gsetdirectories

    maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
    maskpath = [maskroot anim '_' expt];
    %traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
    traceroot = 'i:\Beta\cellTraces\monkey\noZ_control\';
    tracepath = [traceroot anim '_' expt '_cellS'];

    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    if strcmp(anim,'ab2') & strcmp(expt,'u000_087')
        set(G_RChandles.dropTrials,'string','[14:30]')
    else
        set(G_RChandles.dropTrials,'string','[]')
    end

    [maporidiff(ex) mapslopemag(ex)] = Gkernelplots4
    
    %%%%%%%%%
    xshift = 6; yshift = 7;
    prcoverlap = getcontroloverlap(xshift,yshift);
    id = find(prcoverlap > .01);
    TC.kernE{1}(id) = NaN;
    %%%%%%%%%
    
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

SNRcontrol = kernE;
%%

id = find(isnan(SNRcontrol.*SNRactual));
SNRactual(id) = [];
SNRcontrol(id) = [];

hdom = 0:.01:.15;
figure,subplot(3,2,1)
[h hdom] = hist(SNRactual,hdom);  cumh = cumsum(h/sum(h));
bar(hdom,h/length(SNRactual)), xlim([hdom(1) hdom(end)])
title(['med = ' num2str(median(SNRactual)) ';  mu = ' num2str(mean(SNRactual))]), set(gca,'tickdir','out')
subplot(3,2,3)
[h hdom] = hist(SNRcontrol,hdom);  cumhC = cumsum(h/sum(h));
bar(hdom,h/length(SNRcontrol)), xlim([hdom(1) hdom(end)])
title(['med = ' num2str(median(SNRcontrol)) ';  mu = ' num2str(mean(SNRcontrol))]), set(gca,'tickdir','out')

subplot(3,2,5)
plot(hdom,cumh,'k')
hold on
plot(hdom,cumhC)
ylim([0 1.2]),  xlim([hdom(1) hdom(end)])
legend('Actual','Control')

subplot(3,2,2)
hdom = -.5:.04:.5;
Sratio = log10(SNRactual./SNRcontrol);
[h hdom] = hist(Sratio,hdom);
bar(hdom,h/sum(h)), xlim([hdom(1) hdom(end)])
title(['med = ' num2str(median(Sratio)) ';  geomean = ' num2str(10.^mean(Sratio))])
 set(gca,'tickdir','out')

subplot(3,2,4), scatter(SNRactual,SNRcontrol,'.k')
xlabel('SNR (actual)'), ylabel('SNR (control)')
hold on, plot([0 .1],[0 .1])
 
[h p] = ttest2(SNRactual,SNRcontrol)
[h p] = ttest(Sratio)