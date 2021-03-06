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
exptS{2} = 'u009_017'; %Zeus V1; 52 microns more shallow than 9_9
exptS{3} = 'u009_009'; %Zeus V1


% ID all the randpos experiments
%exptSx{1} = 'u009_010';
exptSx{1} = 'u000_106'; %Keith V1
exptSx{2} = 'u009_015'; %Zeus V1; 52 microns more shallow than 9_10
exptSx{3} = 'u009_010'; %Zeus V1

%idExampAll = [6 22 26; 16 20 36; 12 28 57] ;
idExampAll = [6 12 22 26 ;9 29 46 38;5 12 28 57] ;



exid = 2;
%ex = exdom(exid);

exdom = 1:3;
%exdom = 3
clear TC_op PW_op TC_ro PW_ro
%for eid = 1:length(exdom)
for eid = 1:length(exdom)
    
    idExamp = idExampAll(exdom(eid),:);
    
    %idExamp = 1:34
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




%% Accumulate TC info across experiments

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


% Accumlate pairwise info across experiments

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

[TCopAll.D TCopAll.Dnorm TCopAll.rOverlap] = OnOffOverlap(TC_op);


%% Set up layers of the model, starting with the SF domain 

global alpha
alpha = pi; %Number of cycles of SF in sigma*alpha of RF envelope. Used to use 4.
AR = 2; %aspect ratio of input neurons
%clusteringFactor = 2;

%%%%%%%Reset and use low pass cells%%%
%If resetting, make sure this is how maps and everything else is created
resetFlag = 1;
if resetFlag
    TCroAll.sfpref =  TCroAll.sfprefAllpass;
    TCroAll.orisig =  TCroAll.orisigAllpass;
    TCroAll.sfBWLin = TCroAll.fhi - TCroAll.flo;
    TCroAll.sfBW = log2((TCroAll.sfBWLin/2+TCroAll.sfprefAllpass)./TCroAll.sfprefAllpass)*2;
    
    LPthresh = 0.5;
    id = find(~isnan(TCroAll.sfpref) & TCroAll.sfpref>=LPthresh );
    idLP = find(~isnan(TCroAll.sfpref) & TCroAll.sfpref<LPthresh );
else
    id = find(~isnan(TCroAll.sfpref));
    idLP = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize.*TCroAll.flo) );

%id = find(~isnan(TCroAll.sfpref));
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize));
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize) & TCroAll.LPness< 1 );

%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize) & TCroAll.sfpref>0.25 );
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize.*TCroAll.flo));

%TCroAll.sfBWLin = TCroAll.fhi - TCroAll.flo;
%TCroAll.sfBW = log2((TCroAll.sfBWLin/2+TCroAll.sfprefAllpass)./TCroAll.sfprefAllpass)*2;  

%LPthresh = 0.5;
%id = find(~isnan(TCroAll.sfpref) & TCroAll.sfpref>=LPthresh );
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize) & TCroAll.sfpref>=LPthresh );

%idLP = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize) & TCroAll.sfpref<LPthresh );
%idLP = find(~isnan(TCroAll.sfpref) & TCroAll.sfpref<LPthresh );
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize.*TCroAll.flo) );


%%%%%%%%%%%Reset Bandwidth stuff cuz I don't know what I was doing before%%
%The code set some flo values to NaN if they did not meet the 0.61
%threshold for bandpass.  Set these to the lowest sf.

%



%Q.meas.BWlin = (TCroAll.fhi(id) - TCroAll.fhi(id));
%Q.meas.BWlog = (log2(TCroAll.fhi(id)) - log2(flodum(id)))/2;
%Q.meas.BWlog = (log2(TCroAll.fhi(id)) - log2(TCroAll.sfpref(id)));
%Q.meas.BWlog = log2((Q.meas.BWlin+Q.meas.sfpref)./Q.meas.sfpref);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q.sfprefdom = logspace(log10(1/16),log10(10),50)';

%Measured
Q.meas.sfpref = TCroAll.sfpref(id);
Q.meas.size = TCopAll.profileSize(id)/2;

Q.meas.fhi = TCroAll.fhi(id);
Q.meas.flo = TCroAll.flo(id);

for i = 1:length(TCopAll.rawProfile)
    amp(i) = max(TCopAll.rawProfile{i});
end
Q.meas.barAmp = amp(id);
Q.meas.barAmp(find(isnan(Q.meas.size))) = NaN;

%Q.meas.BWlog = log2((Q.meas.sfpref+TCroAll.sfBWLin(id))./Q.meas.sfpref);
Q.meas.BWlin = TCroAll.sfBWLin(id)/2;
Q.meas.BWlog = TCroAll.sfBW(id)/2;
Q.meas.oriBW = TCroAll.orisig(id);
Q.meas.F1F0 = TCroAll.F1F0(id);
Q.meas.Dnorm = TCopAll.Dnorm(id);

%Scale invariance curves; Needs to be done in two domains:
%log domain of SF preference 
Q.SI.size = 1./(alpha*Q.sfprefdom);  %1sig
Q.SI.BWlin = 1./(Q.SI.size*2*pi); %1sig
Q.SI.BWlog = log2(1+alpha/(2*pi)) * ones(size(Q.sfprefdom)); %1sig
Q.SI.oriBW = atan(Q.SI.BWlin./Q.sfprefdom/AR)*180/pi;
%domain of measured SF preference
Q.SI_D.size = 1./(alpha*Q.meas.sfpref);  %1sig
Q.SI_D.BWlin = 1./(Q.SI_D.size*2*pi); %1sig
Q.SI_D.BWlog = log2(1+alpha/(2*pi)) * ones(size(Q.meas.sfpref)); %1sig
Q.SI_D.oriBW = atan(Q.SI_D.BWlin./Q.meas.sfpref/AR)*180/pi;



%%
%Get model parameters by comparing measured data to the SI and DR
%predictions

%compute smearing coefficients
Q.SIsmearXsig = sqrt(nanmedian(Q.meas.size.^2 - Q.SI_D.size.^2)); %pooling window (deg), using scale invariance as input
Q.SIsmearSFsig = sqrt(median(Q.meas.BWlin.^2 - Q.SI_D.BWlin.^2));  %pooling window (c/deg), using scale invariance as input
Q.SIsmearlogSFsig = sqrt(median(Q.meas.BWlog.^2 - Q.SI_D.BWlog.^2)); %pooling window (octaves), using scale invariance as input

%This needs SIsmearXsig, but needs to be done before F1F0 model below
clusteringFactor = getClusteringFactor(Q);  %phase clustering. Set to 1 to make relative phase constant, and absolute phase coupled to the retinotopy

%Model output in artificial SF domain
Q.SI.size_2Cmplx = sqrt(Q.SI.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
Q.SI.BWlin_2Cmplx = sqrt(Q.SI.BWlin.^2 + Q.SIsmearSFsig.^2);  %Model fit of linear SF bw, using scale invariance as inputs
Q.SI.BWlog_2Cmplx = log2((Q.SI.BWlin_2Cmplx+Q.sfprefdom) ./ (Q.sfprefdom));  %Model fit of log SF bw, using scale invariance as inputs
%Q.SI.BWlog_2Cmplx = log2((Q.SI.BWlin_2Cmplx+Q.sfprefdom) ./ (Q.sfprefdom-Q.SI.BWlin_2Cmplx))/2;  %Model fit of log SF bw, using scale invariance as inputs
Q.SI.oriBW_2Cmplx = atan(Q.SI.BWlin_2Cmplx./Q.sfprefdom/AR)*180/pi;
Q.SI.F1F0_2Cmplx = pi*exp(-(Q.SIsmearXsig/clusteringFactor*2*pi*Q.sfprefdom(:)).^2/2); 
Q.SI.F1F0_noCluster_2Cmplx = pi*exp(-(Q.SIsmearXsig*2*pi*Q.sfprefdom(:)).^2/2); 
Q.SI.Dnorm_2Cmplx = zeros(size(Q.sfprefdom(:)));

%Q.SI.BWlog_2Cmplx = sqrt(Q.SI.BWlog.^2 + Q.SIsmearlogSFsig.^2);  %Model fit of log SF bw, using scale invariance as inputs

%Model output in data SF domain
Q.SI_D.size_2Cmplx = sqrt(Q.SI_D.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
Q.SI_D.BWlin_2Cmplx = sqrt(Q.SI_D.BWlin.^2 + Q.SIsmearSFsig.^2);  %Model fit of linear SF bw, using scale invariance as inputs
Q.SI_D.BWlog_2Cmplx = log2((Q.SI_D.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref));  %Model fit of log SF bw, using scale invariance as inputs
%Q.SI_D.BWlog_2Cmplx = log2((Q.SI_D.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref-Q.SI_D.BWlin_2Cmplx))/2;  %Model fit of log SF bw, using scale invariance as inputs
Q.SI_D.oriBW_2Cmplx = atan(Q.SI_D.BWlin_2Cmplx./Q.meas.sfpref/AR)*180/pi;
Q.SI_D.F1F0_2Cmplx = pi*exp(-(Q.SIsmearXsig/clusteringFactor*2*pi*Q.meas.sfpref(:)).^2/2); %analytical solution
Q.SI_D.F1F0_noCluster_2Cmplx = pi*exp(-(Q.SIsmearXsig*2*pi*Q.meas.sfpref(:)).^2/2); %analytical solution
Q.SI_D.Dnorm_2Cmplx = zeros(size(Q.meas.sfpref(:)));

Q
%% Marginal Stats

statout = getAllStats(Q);

%F statistic
[h p] = vartest2(log2(Q.meas.size),log2(Q.meas.sfpref))
[h p] = vartest2(log2(Q.meas.BWlin),log2(Q.meas.sfpref))

%comparisons of variance accounted for between pooling model and linear
%equation
a = log2(Q.meas.sfpref);
b = log2(Q.meas.size);
%% Check fits generated in 'getstats'
linmodwidth =  statout.invwidth.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.invwidth.lineFit_sf_loglog(2);
linmodwidth = 2.^(-linmodwidth);
poolmodwidth = Q.SI_D.size_2Cmplx;

figure
subplot(1,6,1),loglog(Q.meas.sfpref,linmodwidth,'.k'), hold on, loglog(Q.meas.sfpref,poolmodwidth,'.r')
hold on, loglog(Q.meas.sfpref,Q.meas.size,'o')

a = log2(linmodwidth); b = log2(Q.meas.size);
1-var(a-b)/var(b)
[r p] = corrcoef(a,b)
a = log2(poolmodwidth); b = log2(Q.meas.size);
1-var(a-b)/var(b)
[r p] = corrcoef(a,b)

%%%%%linear bw%%%%%%%%%
linmodwidth =  statout.BWlin.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.BWlin.lineFit_sf_loglog(2);
linmodwidth = 2.^(linmodwidth);
poolmodwidth = Q.SI_D.BWlin_2Cmplx;

subplot(1,6,2),loglog(Q.meas.sfpref,linmodwidth,'.k'), hold on, loglog(Q.meas.sfpref,poolmodwidth,'.r')
hold on, loglog(Q.meas.sfpref,Q.meas.BWlin,'o')

a = log2(linmodwidth); b = log2(Q.meas.BWlin);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
a = log2(poolmodwidth); b = log2(Q.meas.BWlin);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)

%%%%%log bw%%%%%%%%%%%%%
logmodwidth =  statout.BWlog.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.BWlog.lineFit_sf_loglog(2);
poolmodwidth = Q.SI_D.BWlog_2Cmplx;

subplot(1,6,3),semilogx(Q.meas.sfpref,logmodwidth,'.k'), hold on, semilogx(Q.meas.sfpref,poolmodwidth,'.r')
hold on, semilogx(Q.meas.sfpref,Q.meas.BWlog,'o')

a = log2(linmodwidth); b = log2(Q.meas.BWlin);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
a = log2(poolmodwidth); b = log2(Q.meas.BWlin);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)


%%%%%ori bw%%%%%%%%%%%%%
orimodwidth =  statout.oriBW.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.oriBW.lineFit_sf_loglog(2);
orimodwidth = 2.^(orimodwidth);
poolmodwidth = Q.SI_D.oriBW_2Cmplx;

subplot(1,6,4),loglog(Q.meas.sfpref,orimodwidth,'.k'), hold on, loglog(Q.meas.sfpref,poolmodwidth,'.r')
hold on, loglog(Q.meas.sfpref,Q.meas.oriBW,'o')

a = log2(orimodwidth); b = log2(Q.meas.oriBW);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
a = log2(orimodwidth); b = log2(Q.meas.oriBW);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)


%%%%%F1F0%%%%%%%%%%%%%
logmodwidth =  statout.F1F0.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.F1F0.lineFit_sf_loglog(2);
poolmodwidth = Q.SI_D.F1F0_2Cmplx;

subplot(1,6,5),semilogx(Q.meas.sfpref,logmodwidth,'.k'), hold on, semilogx(Q.meas.sfpref,poolmodwidth,'.r')
hold on, semilogx(Q.meas.sfpref,Q.meas.F1F0,'o')

a = log2(linmodwidth); b = log2(Q.meas.F1F0);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
a = log2(poolmodwidth); b = log2(Q.meas.F1F0);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)

%%%%%Dnorm%%%%%%%%%%%%%
logmodwidth =  statout.Dnorm.lineFit_sf_loglog(1).*log2(Q.meas.sfpref) + statout.Dnorm.lineFit_sf_loglog(2);
poolmodwidth = Q.SI_D.Dnorm_2Cmplx;

subplot(1,6,6),semilogx(Q.meas.sfpref,logmodwidth,'.k'), hold on, semilogx(Q.meas.sfpref,poolmodwidth,'.r')
hold on, semilogx(Q.meas.sfpref,Q.meas.Dnorm,'o')

a = log2(linmodwidth); b = log2(Q.meas.Dnorm);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
a = log2(poolmodwidth); b = log2(Q.meas.Dnorm);
[r p] = corrcoef(a,b)
1-var(b-a)/var(b)
%%
%%%%%Log Bandwidth%%%%%%%%%%%
plotdom = [1/2 1 2 4];
figure
semilogx(Q.meas.sfpref,Q.meas.BWlog,'.k');
hold on
semilogx(Q.sfprefdom,Q.SI.BWlog,'k');

semilogx(Q.sfprefdom,Q.SI.BWlog_2Cmplx,'g')
set(gca,'XTick',[1/4 1/2 1 2 4 8])
%set(gca,'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{''; '1/2'; '1'; '2'; '4'; '8'})

%hold on
%loglog(sfpref_sim,BWVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)

xlabel('sfreq')
ylabel('log BW')

ylim([1/32 3]), xlim([LPsf 8])
axis square
%xlim(([plotdom(1) plotdom(end)]))


%% Plot

plotdom = [1/32 1/16 1/8 1/4 1/2 1];
LPsf = 1/3;
plotSizeVPrediction2(Q.meas.size,Q.SI_D.size,plotdom); %Measured size vs. scale invariant prediction

subplot(2,2,1), cla
loglog(Q.meas.sfpref,Q.meas.size,'.k')

subplot(2,2,1)
hold on
loglog(Q.sfprefdom,Q.SI.size,'k')
hold on
loglog(Q.sfprefdom,Q.SI.size_2Cmplx,'g')
hold on
loglog(ones(1,length(idLP))*LPsf,TCopAll.profileSize(idLP)/2,'.k') %plot the low pass cells

% id = find(~isnan(TCroAll.sfpref));
% plotSizeVPrediction2(Rsize(id),getSizefromSFpref(sf2sigModel,TCroAll.sfpref(id)),dom)

subplot(2,2,1)
hold on,
%loglog(Q.SI.size,Q.DR.size,'--r','LineWidth',2);
%loglog(.25./sfpref_sim,XsigVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)
set(gca,'XTick',[1/4 1/2 1 2 4 8],'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{''; '1/2'; '1'; '2'; '4'; '8'})
set(gca,'YTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
ylabel('RF width; deg (sigma)')
xlabel('SF preference (c/deg)')

ylim([1/32 1]), xlim([LPsf 8])
axis square

subplot(2,2,2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('RF width; deg (sigma)')
ylabel('# neurons')
xlim(log2([plotdom(1) plotdom(end)]))

subplot(2,2,3)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted width; deg (sigma) 1/(4*sfreq)')
ylabel('RF width; deg (sigma)')
xlim(log2([plotdom(1) plotdom(end)]))


%%%%%Linear Bandwidth%%%%%%%%%%%
plotdom = [1/4 1/2 1 2 4 8];

plotSizeVPrediction2(Q.meas.BWlin,Q.SI_D.BWlin,plotdom)
subplot(2,2,1), cla

subplot(2,2,1), cla
loglog(Q.meas.sfpref,Q.meas.BWlin,'.k')

subplot(2,2,1)
hold on
loglog(Q.sfprefdom,Q.SI.BWlin,'k')
hold on
loglog(Q.sfprefdom,Q.SI.BWlin_2Cmplx,'g')
hold on
loglog(Q.sfprefdom,1./(2*pi*Q.SI.size_2Cmplx))
hold on
loglog(ones(1,length(idLP))*LPsf,TCroAll.sfBWLin(idLP)/2,'.k')
xlim([LPsf 8])

%loglog(2*Q.sfprefdom/pi,Q.DR.BWlin,'--r','LineWidth',2);
set(gca,'XTick',plotdom,'YTick',plotdom)
%loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
xlabel('SF preference (c/deg)'), ylabel('bandwidth; c/deg (sigma)') 

subplot(2,2,2)
%set(gca,'XTickLabel',{'1/4'; '1/2'; '1'; '2'; '4'; '8'})
xlabel('bandwidth; c/deg (sigma)')
ylabel('# neurons')
xlim(log2([plotdom(1) plotdom(end)]))

subplot(2,2,3)
%set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted Bandwidth; c/deg (sigma): sfreq/(alpha*pi*2)')
ylabel('# neurons')
xlim(log2([plotdom(1) plotdom(end)]))


% smearsig = sqrt(mean(log2(1+Q.meas.BWlin./Q.meas.sfpref).^2 - log2(1+2/pi).^2))
% smearsig_cycles = (2.^smearsig-1).*Q.sfprefdom
% smearsig = smearsig.*ones(size(Q.sfprefdom));
% bwE = (2.^sqrt( log2(1+2/pi).^2 + smearsig.^2 ) - 1) .* Q.sfprefdom;
% subplot(2,2,1)
% hold on 
% loglog(2*Q.sfprefdom/pi,bwE)

%%%%%Log Bandwidth%%%%%%%%%%%

%Fit line in log BW domain, then convert back and put it in these axes 


[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),Q.meas.BWlog,[-1 1]);

sigE = (2.^(xhat(1)*log2(Q.sfprefdom)+xhat(2)) - 1)  .*Q.sfprefdom; %algebraic manipulation of Eq. A above gives BW(sfpref)
subplot(2,2,1)
hold on
plot(Q.sfprefdom,sigE,'k')


% y = log2(Q.meas.BWlin./Q.meas.sfpref + 1);
% H = [log2(Q.meas.sfpref) ones(size(Q.meas.sfpref))];
%xhat = inv(H'*H)*H'*y;

% figure,plot(H(:,1),y,'.')
% hold on
% plot(Hhat(:,1),yhat)

[r p] = corrcoef(log2(Q.meas.sfpref),Q.meas.BWlog);
Hhat = [log2(Q.sfprefdom) ones(size(Q.sfprefdom))];
yhat =  Hhat*xhat(:);
figure,plot(log2(Q.meas.sfpref),Q.meas.BWlog,'.k')
hold on
plot(Hhat(:,1),yhat)
hold on
plot(log2(Q.sfprefdom),Q.SI.BWlog_2Cmplx,'g')

ylabel('logBW')
xlabel('logSF')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])

%%%%%Orientation Bandwidth%%%%%%%%%%%


figure
loglog(Q.meas.sfpref,Q.meas.oriBW,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
hold on
loglog(ones(1,length(idLP))*LPsf,TCroAll.orisig(idLP),'.k')

[r p] = corrcoef(log2(Q.meas.sfpref),log2(Q.meas.oriBW));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([LPsf 8.5])
ylim([8 70])
axis square

[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),log2(Q.meas.oriBW),[-.5 1]);
hold on

loglog(Q.sfprefdom,2.^(log2(Q.sfprefdom)*xhat(1) + xhat(2)));
title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(xhat(1)) ])



BWphaseAligned = 1./(2*pi*Q.SI.size_2Cmplx);
oripred_phaseAligned = atan(BWphaseAligned./Q.sfprefdom/AR)*180/pi;

hold on,
loglog(Q.sfprefdom,Q.SI.oriBW_2Cmplx,'g')
hold on,
loglog(Q.sfprefdom,Q.SI.oriBW,'k')
hold on,
loglog(Q.sfprefdom,oripred_phaseAligned,'b')
hold on
loglog(ones(1,length(idLP))*LPsf,TCroAll.orisig(idLP),'.k')

% [sig_x sig_ori] = SFxyPoolingModel(Q.sfprefdom,3,Q.SIsmearSFsig,AR,1100);
% hold on
% loglog(Q.sfprefdom,sig_ori(1,:),'k') 


%% Phase modulation stuff

%[F1F0model] = getF1F0vsSFmodel(Q.sfprefdom(:),Q.SIsmearXsig);
%clusteringFactor = 2;
%F1F0model = pi*exp(-(Q.SIsmearXsig/clusteringFactor*2*pi*Q.sfprefdom(:)).^2/2); %analytical solution

co = 1+eps;
%co = 0.61;

idBP = find(TCroAll.LPness<=co & ~isnan(TCroAll.sfpref));

%idBP = find(TCroAll.LPness<=co & ~isnan(TCroAll.sfpref.*TCroAll.F1F0.*TCroAll.sfBWLin) & ~isnan(TCopAll.Dnorm(:))); %get bandpass population

DnormQuartiles = [prctile(TCopAll.Dnorm(idBP),0) prctile(TCopAll.Dnorm(idBP),25) prctile(TCopAll.Dnorm(idBP),50) prctile(TCopAll.Dnorm(idBP),75) prctile(TCopAll.Dnorm(idBP),100)]
F1F0Quartiles = [prctile(TCroAll.F1F0(idBP),0) prctile(TCroAll.F1F0(idBP),25) prctile(TCroAll.F1F0(idBP),50) prctile(TCroAll.F1F0(idBP),75) prctile(TCroAll.F1F0(idBP),100)]

[r p] = corrcoef(TCroAll.F1F0(idBP),TCopAll.Dnorm(idBP))

figure,
subplot(2,2,1)
semilogx(TCroAll.sfpref(idBP),TCroAll.F1F0(idBP),'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
xlim([.25 8])
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency preference')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.sfpref(idBP),TCroAll.F1F0(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.sfprefdom,Q.SI.F1F0_2Cmplx,'g--')
hold on, semilogx(Q.sfprefdom,Q.SI.F1F0_noCluster_2Cmplx,'g')
xlim([.25 8])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})

subplot(2,2,2)
semilogx(TCroAll.sfBWLin(idBP)/2,TCroAll.F1F0(idBP),'.k')
xlim([.5 4])
%hold on, semilogx(BWLinVsf_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency BW (sigma)')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.F1F0(idBP),TCroAll.sfBWLin(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.SI.BWlin_2Cmplx,Q.SI.F1F0_2Cmplx,'--g')
hold on, semilogx(Q.SI.BWlin_2Cmplx,Q.SI.F1F0_noCluster_2Cmplx,'g')
%hold on, semilogx(Q.SI.BWlin,F1F0model,'g')
%xlim([.6 3])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})


subplot(2,2,3)
id = find(~isnan(TCopAll.profileSize(idBP)));
semilogx(TCroAll.sfpref(idBP(id)).*TCopAll.profileSize(idBP(id))*pi/2,TCroAll.F1F0(idBP(id)),'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
xlim([.25 8])
set(gca,'XTick',[1 2 4 8 16])
xlabel('measured RF width / SI prediction')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.sfpref(idBP(id)).*TCopAll.profileSize(idBP(id))*pi,TCroAll.F1F0(idBP(id)));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.sfprefdom.*Q.SI.size_2Cmplx*pi,Q.SI.F1F0_2Cmplx,'--g')
hold on, semilogx(Q.sfprefdom.*Q.SI.size_2Cmplx*pi,Q.SI.F1F0_noCluster_2Cmplx,'g')
xlim([.5 16])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})

subplot(2,2,4)
semilogx(TCroAll.sfBWLin(idBP)/2./TCroAll.sfpref(idBP)*2,TCroAll.F1F0(idBP),'.k')
xlim([.5 8])
%hold on, semilogx(BWLinVsf_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('measured bandwidth / SI prediction')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.F1F0(idBP),TCroAll.sfBWLin(idBP)/2./TCroAll.sfpref(idBP)*2);
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.SI.BWlin_2Cmplx./Q.sfprefdom*2,Q.SI.F1F0_2Cmplx,'--g')
hold on, semilogx(Q.SI.BWlin_2Cmplx./Q.sfprefdom*2,Q.SI.F1F0_noCluster_2Cmplx,'g')
%hold on, semilogx(Q.SI.BWlin,Q.SI.F1F0_2Cmplx,'g')
%xlim([.6 3])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})


%%
%compute variance accounted for
f = log2(Q.meas.oriBW);

ffit = log2(oripred_D);  %Use the pooling model as model
varacc_pool = (var(f(:))-var(f(:)-ffit(:)))/var(f(:))
mse_pool = sqrt(mean((f-ffit).^2))
ffit = xhat(1)*log2(Q.meas.sfpref) + xhat(2); %Use linear fit for the model
varacc_lin = (var(f(:))-var(f(:)-ffit(:)))/var(f(:))
mse_lin = sqrt(mean((f-ffit).^2))
mse_pool/mse_lin

subplot(3,1,2)
loglog(Q.meas.BWlin,Q.meas.oriBW,'.k'), xlabel('linear SF bandwidth  (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(Q.meas.BWlin),log2(Q.meas.oriBW));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 8.5])
ylim([8 45])

[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.BWlin),log2(Q.meas.oriBW),[-.5 1]);
hold on
loglog(2.^dom,2.^ffit)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(xhat(1)) ])


oripred = atan(Q.SI.BWlin./Q.sfprefdom/AR)*180/pi;
%oripred = atan(Q.SI.BWlin_2Cmplx./Q.sfprefdom/AR)*180/pi;

hold on,
plot(Q.SI.BWlin,oripred)


subplot(3,1,3)
semilogy(Q.meas.BWlog,Q.meas.oriBW,'.k'), xlabel('logarithmic SF bandwidth  (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef((Q.meas.BWlog),log2(Q.meas.oriBW));
set(gca,'XTick',[linspace(.25,3,6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 3])
ylim([8 45])

[xhat ffit dom varaccount] = LinetotalLS((Q.meas.BWlog),log2(Q.meas.oriBW),[2 2]);
hold on
semilogy(dom,2.^ffit)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(xhat(1)) ])


oripred = atan(Q.SI.BWlin./Q.sfprefdom/AR)*180/pi;
%oripred = atan(Q.SI.BWlin_2Cmplx./Q.sfprefdom/AR)*180/pi;

hold on,
plot(Q.SI.BWlog,oripred)



%% Plot of orientation selectivity vs spatial freq preference
global DM
clear kernrawAll
id1 = 1;
for i = 1:length(TC_ro)
   
    idx = id1:id1+size(TC_ro{i}.tcorisfraw{1},1)-1;
    kernrawAll(idx,:,:) = TC_ro{i}.tcorisfraw{1};
    
    id1 = idx(end)+1;
    
end

plotOrivSF(kernrawAll,TCroAll.opref,TCroAll.sfpref,DM.oridom,DM.sfdom,[.5 1 2 8])



%% Pairwise stuff
Dmax = 525;
Dbins = 0:75:Dmax;

e = 2;
SFOriPhaseClustering(PW_ro{2},Dbins)

sizeClustering(PW_op{2},Dbins)

sizeClustering(PWopAll,Dbins)
SFOriPhaseClustering(PWroAll,Dbins)

RetvPhaseClustering(PWroAll,PWopAll)

%Compute relative spatial phase by subtracting the retinotopy
musf = sqrt(PWroAll.dsfpair(:,1).*PWroAll.dsfpair(:,2));
dcycles = abs(PWroAll.dphase)/360;
ddegrees = dcycles./musf;
PWroAll.dcycles = ddegrees - PWopAll.dpos;


%% Assessing response amplitude to bars vs. SF tuning

nanmedian(Q.meas.flo)
nanmedian(Q.meas.fhi)

idLow = (find(Q.meas.fhi<2.5));
idMid = (find(Q.meas.fhi>2.5 & Q.meas.flo<2.5));
idHi = (find(Q.meas.flo>2.5));

Ncells = length(~isnan(Q.meas.fhi));

length(idLow)/Ncells
length(idMid)/Ncells
length(idHi)/Ncells

% LowMedian = nanmedian(TCopAll.amp(idLow))
% HighMedian = nanmedian(TCopAll.amp(idnLow))
% HighMedian/LowMedian
% [h p] = ttest2(TCopAll.amp(idLow),TCopAll.amp(idnLow))

LowMedian = nanmedian(Q.meas.barAmp(idLow))
MidMedian = nanmedian(Q.meas.barAmp(idMid))
HighMedian = nanmean(Q.meas.barAmp(idHi))
PercAttenuation_Hi = HighMedian/MidMedian
PercAttenuation_Lo = LowMedian/MidMedian

[h p] = ttest2(Q.meas.barAmp(idLow),Q.meas.barAmp(idMid))
[h p] = ttest2(Q.meas.barAmp(idHi),Q.meas.barAmp(idMid))

id = find(~isnan(Q.meas.sfpref(:).*Q.meas.barAmp(:)));
[r p] = corrcoef(Q.meas.sfpref(id),Q.meas.barAmp(id))



%% Get simulation of complex cell generation

global alpha
alpha = pi;

%sf2sigModel = 'NL sat';
%sf2sigModel = 'exp decay';
sf2sigModel = 'classic'

sfBands = linspace(1,5,5);
if strcmp(sf2sigModel,'NL sat') 
    xsigDom = .2;
    sfsigDom = 1;
elseif strcmp(sf2sigModel,'classic') 
    xsigDom = Q.SIsmearXsig;
    sfsigDom = Q.SIsmearSFsig;
end
%sfsigDom = .5
%[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation_map(sfBands,xsigDom,sfsigDom,sf2sigModel);

[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation(sfBands,xsigDom,sfsigDom,sf2sigModel);
%[XsigVsf_sim BWLinVsf_sim sfpref_sim F1vsBand_sim orisigVsf] = V1_ComplexCell_generation_2Dx(sfBands,xsigDom,sfsigDom,sf2sigModel);

BWLinVsf_sim = BWLinVsf_sim*2; %Make it 2sigma to be consistant with stuff below.
%sfpref_sim = sfBands

%% DATA yield stats

N = length(TCroAll.sfpref);

%random grating yield
idRGbad = find(isnan(TCroAll.sfprefAllpass)); %w/o low pass exclusion
idRGbadLP = find(isnan(TCroAll.sfpref)); %w low pass exclusion
RGexclusion = length(idRGbad)/N  %Exclusion from fitting alone
NLP = length(idRGbadLP)-length(idRGbad)  %No of Low pass cells that passed fitting criterion

%random bar yield
idRObad = find(isnan(TCopAll.profileSize));
ROexclusion = length(idRObad)/N

idAllbad = unique([idRObad; idRGbadLP]);
Allexclusion = length(idAllbad)/N

