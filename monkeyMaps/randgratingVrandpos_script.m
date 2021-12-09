pF0
pRev

global G_RChandles G_handles TC PW idExamp

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
%%
clear animS exptS exptSx

animS{1} = 'ax2';
animS{2} = 'ax3';
animS{3} = 'ax3';

% ID all the random orisf experiments
%exptS{1} = 'u009_009';
exptS{1} = 'u000_108'; %Keith V1
exptS{2} = 'u009_017'; %Zeus V1
exptS{3} = 'u009_009'; %Zeus V1


% ID all the randpos experiments
%exptSx{1} = 'u009_010';
exptSx{1} = 'u000_106'; %Keith V1
exptSx{2} = 'u009_015'; %Zeus V1
exptSx{3} = 'u009_010'; %Zeus V1

%idExampAll = [6 22 26; 16 20 36; 12 28 57] ;
idExampAll = [6 12 22 26 ;5 16 20 36 ;5 12 28 57] ;



%exid = 2;
%ex = exdom(exid);

exdom = 1:3;
%exdom = 3
clear TC_op PW_op TC_ro PW_ro
%for eid = 1:length(exdom)
    
%exdom = 2;
for eid = 1:length(exdom)
    
    idExamp = idExampAll(exdom(eid),:);
    
    ex = exdom(eid);
    
    anim = animS{ex};
    expt = exptS{ex};
    exptx = exptSx{ex};
    
    %%%%%load randori mask%%%%%%%
    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    load(maskpath,'maskS')
    
    %%%%%load oripos traces and process%%%%%%%
    traceroot = 'C:\CellTraces\randpos_alignedto_randori\';
    tracepath = [traceroot anim '_' exptx(1:8) '_cellS'];
    load(tracepath,'cellS')
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' exptx(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' exptx(1:8)])
    Gsetdirectories
    set(G_RChandles.dropTrials,'string','[]')
    Grandposplots3
    TC_op{eid} = TC;
    PW_op{eid} = PW;
    
    
    %%%%%%load randori traces and process%%%%%%%%
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    load(tracepath,'cellS')
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    set(G_RChandles.dropTrials,'string','[]')
    if strcmp(expt,'u009_009') & strcmp(anim,'ax3')
        set(G_RChandles.dropTrials,'string','[21:65]');
    end
    Gkernelplots4
    
    TC_ro{eid} = TC;
    PW_ro{eid} = PW;
    
end




%% Accumlate TC info across experiments

fnames = fieldnames(TC_ro{1});
for e = 1:length(TC_ro)
    for i = 1:length(fnames)
        if e == 1
            eval(['TCroAll.' fnames{i} ' = [];']);
        end
        
        dum = eval(['TC_ro{' num2str(e) '}.'  fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval(['TCroAll.' fnames{i} '= [TCroAll.' fnames{i} '; dum(:)];']);
            %         else
            %             eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
end


fnames = fieldnames(TC_op{1});
for e = 1:length(TC_op)
    for i = 1:length(fnames)
        if e == 1
            eval(['TCopAll.' fnames{i} ' = [];']);
        end
        
        dum = eval(['TC_op{' num2str(e) '}.'  fnames{i} '{1}{1};']); dim = size(dum);
        if min(dim) == 1
            eval(['TCopAll.' fnames{i} '= [TCopAll.' fnames{i} '; dum(:)];']);
            %         else
            %             eval([fnames{i} '= [' fnames{i} '; dum];']);
        end
    end
end

%animID = [animID ex*ones(1,length(TC.xpos{1}))];  %for bootstrap later


%% Accumlate pairwise info across experiments

fnames = fieldnames(PW_ro{1});
for e = 1:length(PW_ro)
    for i = 1:length(fnames)
        if e == 1
            eval(['PWroAll.' fnames{i} ' = [];']);
        end
        
        dum = eval(['PW_ro{' num2str(e) '}.'  fnames{i} '{1};']); dim = size(dum);
        if min(dim) == 1
            eval(['PWroAll.' fnames{i} '= [PWroAll.' fnames{i} '; dum(:)];']);
        else
            eval(['PWroAll.' fnames{i} '= [PWroAll.' fnames{i} '; dum];']);
        end
    end
end


fnames = fieldnames(PW_op{1});
for e = 1:length(PW_op)
    for i = 1:length(fnames)
        if e == 1
            eval(['PWopAll.' fnames{i} ' = [];']);
        end
        
        dum = eval(['PW_op{' num2str(e) '}.'  fnames{i} '{1}{1};']); dim = size(dum);
        if min(dim) == 1
            eval(['PWopAll.' fnames{i} '= [PWopAll.' fnames{i} '; dum(:)];']);
        else
            eval(['PWopAll.' fnames{i} '= [PWopAll.' fnames{i} '; dum];']);
        end
    end
end

%animID = [animID ex*ones(1,length(PW.xpos{1}))];  %for bootstrap later

%% Get simulation of complex cell generation

sfBands = linspace(.5,6,5);
xsigDom = .22;
sfsigDom = 2;

[XsigVsf BWLinVsf, BWVsf, sfpref] = V1_ComplexCell_generation(sfBands,xsigDom,sfsigDom);


%% Compare size to predicted size from spatial frequency preference

%e = 1;
%sf = TC_ro{e}.sfpref{1};
%sfpred = 1./(2*TC_op{e}.profileSize{1}{1});

Rsize = TCopAll.profileSize; %actual RF width
Rsize = TCopAll.xsize

Rsizepred = 1./(2*TCroAll.sfpref); %predicted width (2sigma) from spatial frequency
%Rsizepred = TCroAll.RFsizeFromIFFT

% param = [1.76 0.57 0.197];
% Rsizepred = (param(1)*exp(-TCroAll.sfpref*param(2)) + param(3))/2;



plotSizeVPrediction(Rsize,Rsizepred)


for i = 1:length(sfsigDom)
    
    cdom = [1 0 0; 0 1 0; 0 0 1];
    for j = 1:length(xsigDom)
        subplot(2,2,1)
        loglog(1./(sfpref(:,j,i))/2,2*XsigVsf(:,j,i),'o-','LineWidth',2,'MarkerSize',5,'Color',cdom(j,:)/i)
        set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
        ylabel('measured RF width')
        xlabel('predicted RF width (2sig): .5/(sfpref)')
        %xlim([.05 1.5]), ylim([.05 1.5])
        hold on
        
        %Plug in fit from Ringach data
        sfprefdom = logspace(log10(1/16),log10(8),30);
        param = [1.76 0.57 0.197];
        sigE = param(1)*exp(-sfprefdom*param(2)) + param(3);
        sigE = sigE/4;  %model was fit using 4sig.
        loglog(.5./sfprefdom,sigE*2,'b');
        
    end
    
    clear leg
    for i = 1:length(xsigDom)
        leg{i} = num2str(xsigDom(i))
    end
    legend(leg)
end

subplot(2,2,3), xlabel('predicted size from sf pref (size = .5/(sfpref)')

%% Compute pooling window in mm of cortical space

param = [1.76 0.57 0.197];
sigE = (param(1)*exp(-TCroAll.sfpref*param(2)) + param(3))/4; %Model from Ringach data
Rsizepred = sigE*2; %2sigma
Rsize = TCopAll.profileSize; %actual RF width
%Rsizepred = 1./(2*TCroAll.sfpref);

id = find(~isnan(Rsize.*Rsizepred));
Rsize = Rsize(id);
Rsizepred = Rsizepred(id);

figure,
subplot(1,2,1)
plot((Rsizepred/2).^2,(Rsize/2).^2,'o');
xlabel('sigPrediction^2'), ylabel('sigMeasured^2')
hold on, plot([0 .3],[0 .3])

[r p] = corrcoef(Rsize.^2,Rsizepred.^2)
title(['r=', num2str(r(1,2)) ' p=' num2str(p(1,2))])

subplot(1,2,2),
id = find(Rsize>Rsizepred);
smearSig = sqrt((Rsize(id)/2).^2 - (Rsizepred(id)/2).^2);
hist(real(smearSig))
xlabel('smear = sqrt(sigMeas^2 - sigPred^2)')
muSmear = 2*mean(smearSig); %degrees 2sigma
title(['Smear (2sig)=' num2str(muSmear)])

magFac = 1; %mm/deg
corticalSmear = muSmear*magFac; %mm 2sigma

%% Compare size to size prediction from BW

%e = 1;
%sfBW = TC_ro{e}.sfBWLin{1}/2; %HPc - LPc (@ 1/sqrt(max)) (~2sig) 
%sfBWpred = 1./(pi*TC_op{e}.profileSize{1}{1}/2); %(BW is 2sig); profile size is 2sig

%sig_x = 1/(2*pi*sig_f) 
%~1/4sfpref... sig_f = 2sfpref/pi

Rsize = TCopAll.profileSize; %actual RF width; 2sigma

%sfprefEst = TCroAll.sfBWLin*pi/4;
%Rsizepred = 1./(2*sfprefEst); %2sig size estimate from BW
Rsizepred = 1./(TCroAll.sfBWLin*pi); %2sig size estimate from 2sig BW
%Rsizepred = 1./((TCroAll.fhi-TCroAll.flo)*pi); %2sig size estimate from 2sig BW
%Rsizepred = 1./((TCroAll.fhi-TCroAll.sfpref)*2*pi); %2sig size estimate from 2sig BW

plotSizeVPrediction(Rsize,Rsizepred)

for i = 1:length(sfsigDom)
    
    cdom = [1 0 0; 0 1 0; 0 0 1];
    for j = 1:length(xsigDom)
        
        subplot(2,2,1)
        %BWLinVsf from simulation is 2sigma.  
        loglog(2./(2*pi*BWLinVsf(:,j,i)),2*XsigVsf(:,j,i),'o-','LineWidth',2,'MarkerSize',5,'Color',cdom(j,:)/i)
        set(gca,'XTick',[1/16 1/8 .25 .5 1 2],'YTick',[1/16 1/8 .25 .5 1 2])
        ylabel('measured RF width')
        xlabel('predicted RF width (2sig): 2/(2*pi*BWsig)')
        %xlim([.05 1.5]), ylim([.05 1.5])
        hold on
        
        %Plug in fit from Ringach data
        sfprefdom = logspace(log10(1/16),log10(8),30);
        param = [1.76 0.57 0.197];
        sigE = param(1)*exp(-sfprefdom*param(2)) + param(3);
        sigE = sigE/4;  %model was fit using 4sig.
        loglog(.5./sfprefdom,sigE*2,'b');
        

    end
    
    clear leg
    for i = 1:length(xsigDom)
        leg{i} = num2str(xsigDom(i));
    end
    legend(leg)
end
subplot(2,2,3), xlabel('predicted size from sf BW (size = 2/(2BWpi)')

%% Compute pooling window in cyc/deg of cortical space. i.e. the smearing function that creates wider sf tuning

param = [1.76 0.57 0.197];
sigE = (param(1)*exp(-TCroAll.sfpref*param(2)) + param(3))/4; %Model from Ringach data
BWpred = 2./(sigE*2*pi); %2sigma
BW = TCroAll.sfBWLin; %actual BW
%BWpred = 4*TCroAll.sfpref/pi; %2sigma from 4sig_x = 1/sfpref

id = find(~isnan(BW.*BWpred));
BW = BW(id);
BWpred = BWpred(id);

figure,
subplot(1,2,1)
plot((BWpred/2).^2,(BW/2).^2,'o');
xlabel('sigPrediction^2'), ylabel('sigMeasured^2')
hold on, plot([0 3],[0 3])

[r p] = corrcoef(BW.^2,BWpred.^2)
title(['r=', num2str(r(1,2)) ' p=' num2str(p(1,2))])

subplot(1,2,2),
id = find(BW>BWpred);
smearSig = sqrt((BW(id)/2).^2 - (BWpred(id)/2).^2);
hist(real(smearSig))
xlabel('smear = sqrt(sigMeas^2 - sigPred^2) c/deg')
muSmear = 2*mean(smearSig); %degrees 2sigma
title(['Smear (2sig)=' num2str(muSmear) 'c/deg'])

magFac = 1; %mm/deg
corticalSmear = muSmear*magFac; %mm 2sigma

%%

figure,
subplot(1,2,1)
scatter(log2(TCroAll.sfpref),log2(TCroAll.sfBW)), axis square, hold on, plot([0 4],[0 4]), xlabel('sfpref'),ylabel('sfhi/sflo')
subplot(1,2,2)
scatter(log2(TCroAll.sfpref),log2(TCroAll.sfBWLin)), axis square, hold on, plot([0 4],[0 4]), xlabel('sfpref'),ylabel('sfhi-sflo')



%%
e = 2;
SFOriPhaseClustering(PW_ro{3})

SFOriPhaseClustering(PWroAll)

RetvPhaseClustering(PWroAll,PWopAll)
%%
figure, 
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));
semilogx(TCroAll.sfpref,TCroAll.orisig,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma)')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),TCroAll.orisig(id));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
xlim([.2 8.5])

H = [log2(TCroAll.sfpref(id)) ones(size(TCroAll.sfpref(id)))];
pa = inv(H'*H)*H'*TCroAll.orisig(id)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'deg/octave'])






