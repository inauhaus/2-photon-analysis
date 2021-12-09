global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp

set(G_handles.datadir,'string','c:\2p_data\')
set(G_handles.analyzedir,'string','c:\2p_data\AnalyzerFiles\')

%set(G_handles.datadir,'string','C:\2p_data\')
%set(G_handles.analyzedir,'string','C:\2p_data\AnalyzerFiles\')

%%%%%%%First, left hemisphere%%%%%%%%%%%

idExamp = [100 217; 211 92; 156 140; 138 191];  %left hemi


%% ab8; u000_097   (~190um deep; shallowest) 

anim = 'ab8';
expt = 'u000_097';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 4;
[oriang{1} orimag{1} sfpref{1} sfmag{1} dpmask{1}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.1');
hh = makeMapFilter;
plotMapExamples(hh)

%% ab8; u000_093   (~250um deep; middle) 

anim = 'ab8';
expt = 'u000_093';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 6.5;
[oriang{2} orimag{2} sfpref{2} sfmag{2} dpmask{2}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.1');
hh = makeMapFilter;
plotMapExamples(hh) 

%% ab8; u000_095   (~310um deep; deepest) 

anim = 'ab8';
expt = 'u000_095';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 4;
[oriang{3} orimag{3} sfpref{3} sfmag{3} dpmask{3}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

%% Compute depth continuity

d1 = 1; d2 = 3;

idROI = find( dpmask{d1}.*dpmask{d2} == 1 );

[oriCont_mu oriCont_sig sfCont_mu sfCont_sig] = depthContinuity(oriang,sfpref,dpmask,[d1 d2])

[rori] = circCorr(oriang{d1}(idROI),oriang{d2}(idROI), 180)
[rsf p] = corrcoef(log2(sfpref{d1}(idROI)),log2(sfpref{d2}(idROI)))

%%%%%%%Now right hemisphere%%%%%%%%%%%


%% ab8; u001_008   (~180um deep; shallower) 
idExamp = [183 167; 223 117; 220 139; 235 132];  %right hemi
clear Tens Tens_var %I didn't save these because they are too big with faster scanning

anim = 'ab8';
expt = 'u001_008';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3;
[oriang{1} orimag{1} sfpref{1} sfmag{1} dpmask{1}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

%% ab8; u001_007   (~230um deep; middle) 

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
[oriang{2} orimag{2} sfpref{2} sfmag{2} dpmask{2}] = OriSf_WideField(dthresh);
plotMapExamples(hh)

%% ab8; u001_009   (~280um deep; deepest) 

anim = 'ab8';
expt = 'u001_009';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3.5;
[oriang{3} orimag{3} sfpref{3} sfmag{3} dpmask{3}] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

%% Compute depth continuity

d1 = 1; d2 = 3;
idROI = find(bw);

[oriCont_mu oriCont_sig sfCont_mu sfCont_sig] = depthContinuity(oriang,sfpref,dpmask,[d1 d2])

[rori] = circCorr(oriang{d1}(idROI),oriang{d2}(idROI), 180)
[rsf p] = corrcoef(log2(sfpref{d1}(idROI)),log2(sfpref{d2}(idROI)))

%%%%%%%%%%%%%Other animals%%%%%%%%%%%%%%


%% ab9; u000_084 (bigger spot, thats decent)

anim = 'ab9';
expt = 'u000_084';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)]) %N.B.  A trial from experiment 80 (aborted) was averaged in with f0m and f0m_var.

dthresh = 1.5;
[oriang orimag sfpref sfmag dpmask] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

% circcorr2D(oriang,dpmask,1)
% circcorr2D(sfpref,dpmask,0)
%% ab9; u000_059 (a nice patch)

anim = 'ab9';
expt = 'u000_059';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)]) 

dthresh = 2;
[oriang orimag sfpref sfmag dpmask] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

% circcorr2D(oriang,dpmask,1)
% circcorr2D(sfpref,dpmask,0)
%% ac0; u000_047  (ok)

anim = 'ac0';
expt = 'u000_047';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1.5;
[oriang orimag sfpref sfmag dpmask] = OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

% circcorr2D(oriang,dpmask,1)
% circcorr2D(sfpref,dpmask,0)
%% ac0; u000_062 (this one is pretty crappy)

anim = 'ac0';
expt = 'u000_062';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['C:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 3;
OriSf_WideField(dthresh);

set(G_handles.Lwidth,'string','.5');
hh = makeMapFilter;
plotMapExamples(hh)

%% Test if pinwheels are biased to Low/high SF

global sfdist SFatPinwheel

SFPinAll = [];
sfdistAll = [];
for i = 1:length(SFatPinwheel)   %loop each ROI
    SFPinAll = [SFPinAll; SFatPinwheel{i}(:)];    
    sfdistAll = [sfdistAll; sfdist{i}(:)];
end

hdom = logspace(log10(.5),log10(8),30);
h1 = hist(SFPinAll,hdom);
h2 = hist(sfdistAll,hdom);

figure,
subplot(2,1,1), bar(h1), ylabel('no of pinwheels')
subplot(2,1,2), bar(h2/sum(h2)), xlabel('percentage of pixels')

[h p]= ttest2(log2(SFPinAll),log2(sfdistAll))