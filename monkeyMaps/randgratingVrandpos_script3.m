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
    
%exdom = 3;
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

sf2sigModel = 'NL sat';
%sf2sigModel = 'exp decay';
sf2sigModel = 'classic'

sfBands = linspace(1,5,5);
if strcmp(sf2sigModel,'NL sat') 
    xsigDom = .2;
    sfsigDom = 1;
elseif strcmp(sf2sigModel,'classic') 
    xsigDom = .23;
    sfsigDom = .78;
end
%sfsigDom = .5
%[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation_map(sfBands,xsigDom,sfsigDom,sf2sigModel);
[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation(sfBands,xsigDom,sfsigDom,sf2sigModel);
%[XsigVsf_sim BWLinVsf_sim sfpref_sim F1vsBand_sim orisigVsf] = V1_ComplexCell_generation_2D(sfBands,xsigDom,sfsigDom,sf2sigModel);
BWLinVsf_sim = BWLinVsf_sim*2; %Make it 2sigma to be consistant with stuff below.
%sfpref_sim = sfBands


%% Set up layers of the model, starting with the SF domain 

id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize));
%id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize.*TCroAll.flo));

Q.sfprefdom = logspace(log10(1/16),log10(10),50)';


%Measured
Q.meas.sfpref = TCroAll.sfpref(id);
Q.meas.size = TCopAll.profileSize(id)/2;
Q.meas.BWlin =  TCroAll.sfBWLin(id)/2;
%Q.meas.BWlin = TCroAll.fhi(id)-TCroAll.sfpref(id);
%Q.meas.BWlog = TCroAll.sfBW(id)/2;
Q.meas.BWlog = log2((Q.meas.sfpref+TCroAll.sfBWLin(id))./Q.meas.sfpref);
Q.meas.BWlog = TCroAll.sfBWsig(id);
Q.meas.oriBW = TCroAll.orisig(id);

%Scale invariance curves; Needs to be done in two domains:
%log domain of SF preference 
Q.SI.size = 1./(4*Q.sfprefdom);  %1sig
Q.SI.BWlin = 1./(Q.SI.size*2*pi); %1sig
Q.SI.BWlog = log2(1+2/pi) * ones(size(Q.sfprefdom)); %1sig
Q.SI.oriBW = atan(Q.SI.BWlin./Q.sfprefdom)*180/pi;
%domain of measured SF preference
Q.SI_D.size = 1./(4*Q.meas.sfpref);  %1sig
Q.SI_D.BWlin = 1./(Q.SI_D.size*2*pi); %1sig
Q.SI_D.BWlog = log2(1+2/pi) * ones(size(Q.meas.sfpref)); %1sig
Q.SI_D.oriBW = atan(Q.SI_D.BWlin./Q.meas.sfpref)*180/pi;

%Simple cell (Ringach data) curves; Needs to be done in two domains:
%log domain of SF preference 
Q.DR.size = getSizefromSFpref(sf2sigModel,Q.sfprefdom);  %1sig
Q.DR.BWlin = 1./(Q.DR.size*2*pi); %1sig. Based on FT of Gaussian
Q.DR.BWlog = log2((Q.sfprefdom+Q.DR.BWlin)./(Q.sfprefdom)); %1sig
Q.DR.oriBW = atan(Q.DR.BWlin./Q.sfprefdom)*180/pi;
%domain of measured SF preference
Q.DR_D.size = getSizefromSFpref(sf2sigModel,Q.meas.sfpref);
Q.DR_D.BWlin = 1./(Q.DR_D.size*2*pi); %1sig
Q.DR_D.BWlog = log2((Q.meas.sfpref+Q.DR_D.BWlin)./(Q.meas.sfpref)); %1sig
Q.DR_D.oriBW = atan(Q.DR_D.BWlin./Q.meas.sfpref)*180/pi;


%%
%Get model parameters by comparing measured data to the SI and DR
%predictions
Q.SIsmearXsig = sqrt(mean(Q.meas.size.^2 - Q.SI_D.size.^2)) %pooling window (deg), using scale invariance as input
Q.DRsmearXsig = sqrt(mean(Q.meas.size.^2 - Q.DR_D.size.^2)) %pooling window (deg), using simple cells as input

Q.SIsmearSFsig = sqrt(mean(Q.meas.BWlin.^2 - Q.SI_D.BWlin.^2))  %pooling window (c/deg), using scale invariance as input
Q.DRsmearSFsig = sqrt(mean(Q.meas.BWlin.^2 - Q.DR_D.BWlin.^2))   %pooling window (c/deg), using simple cells as input

Q.SIsmearlogSFsig = sqrt(mean(Q.meas.BWlog.^2 - Q.SI_D.BWlog.^2)) %pooling window (octaves), using scale invariance as input
Q.DRsmearlogSFsig = sqrt(mean(Q.meas.BWlog.^2 - Q.DR_D.BWlog.^2))  %pooling window (octaves), using simple cells as input

Q.SI.size_2Cmplx = sqrt(Q.SI.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
Q.DR.size_2Cmplx = sqrt(Q.DR.size.^2 + Q.DRsmearXsig.^2); %Model fit of size, using simple cells as inputs

Q.SI.BWlin_2Cmplx = sqrt(Q.SI.BWlin.^2 + Q.SIsmearSFsig.^2);  %Model fit of linear SF bw, using scale invariance as inputs
Q.DR.BWlin_2Cmplx = sqrt(Q.DR.BWlin.^2 + Q.DRsmearSFsig.^2);  %Model fit of linear SF bw, using simple cells as inputs

Q.SI.BWlog_2Cmplx = log2((Q.SI.BWlin_2Cmplx+Q.sfprefdom) ./ (Q.sfprefdom));  %Model fit of log SF bw, using scale invariance as inputs
Q.DR.BWlog_2Cmplx = log2((Q.DR.BWlin_2Cmplx+Q.sfprefdom) ./ (Q.sfprefdom)); %Model fit of log SF bw, using simple cells as input

%Q.SI.BWlog_2Cmplx = sqrt(Q.SI.BWlog.^2 + Q.SIsmearlogSFsig.^2);  %Model fit of log SF bw, using scale invariance as inputs
%Q.DR.BWlog_2Cmplx = sqrt(Q.DR.BWlog.^2 + Q.DRsmearlogSFsig.^2);  %Model fit of log SF bw, using simple cells as input

%%
% 
% Q.SI.size_2Cmplx = sqrt(Q.SI_D.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
% Q.DR.size_2Cmplx = sqrt(Q.DR_D.size.^2 + Q.DRsmearXsig.^2); %Model fit of size, using simple cells as inputs
% 
% Q.SI.BWlin_2Cmplx = 1./(2*pi*Q.SI.size_2Cmplx) ;  %Model fit of linear SF bw, using scale invariance as inputs
% Q.DR.BWlin_2Cmplx = 1./(2*pi*Q.DR.size_2Cmplx);  %Model fit of linear SF bw, using simple cells as inputs
% 
% Q.SIsmearSFsig = sqrt(median(Q.meas.BWlin.^2 - Q.SI.BWlin_2Cmplx.^2))  %pooling window (c/deg), using scale invariance as input
% Q.DRsmearSFsig = sqrt(median(Q.meas.BWlin.^2 - Q.DR.BWlin_2Cmplx.^2))   %pooling window (c/deg), using simple cells as input
% 
% 
% Q.SI.size_2Cmplx = sqrt(Q.SI_D.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
% Q.DR.size_2Cmplx = sqrt(Q.DR_D.size.^2 + Q.DRsmearXsig.^2); %Model fit of size, using simple cells as inputs
% 
% Q.SI.BWlin_2Cmplx = 1./(2*pi*Q.SI.size_2Cmplx) ;  %Model fit of linear SF bw, using scale invariance as inputs
% Q.DR.BWlin_2Cmplx = 1./(2*pi*Q.DR.size_2Cmplx);  %Model fit of linear SF bw, using simple cells as inputs
% 
% %Q.SI.BWlin_2Cmplx = sqrt(Q.SIsmearSFsig.^2 + Q.SI.BWlin_2Cmplx.^2)  %pooling window (c/deg), using scale invariance as input
% %Q.DR.BWlin_2Cmplx = sqrt(Q.DRsmearSFsig.^2 + Q.DR.BWlin_2Cmplx.^2)   %pooling window (c/deg), using simple cells as input
% 
% Q.SI.BWlog_2Cmplx = log2((Q.SI.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref));  %Model fit of log SF bw, using scale invariance as inputs
% Q.DR.BWlog_2Cmplx = log2((Q.DR.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref)); %Model fit of log SF bw, using simple cells as input
% 
% Q.SIsmearlogSFsig = sqrt(median(Q.meas.BWlog.^2 - Q.SI.BWlog_2Cmplx.^2)) %pooling window (octaves), using scale invariance as input
% Q.DRsmearlogSFsig = sqrt(median(Q.meas.BWlog.^2 - Q.DR.BWlog_2Cmplx.^2))  %pooling window (octaves), using simple cells as input






%%
%%%%%Log Bandwidth%%%%%%%%%%%
plotdom = [1/2 1 2 4];
figure
subplot(2,2,1)
semilogx(Q.meas.sfpref,Q.meas.BWlog,'.k');
hold on
semilogx(Q.sfprefdom,Q.DR.BWlog,'--k');
hold on
semilogx(Q.sfprefdom,Q.SI.BWlog,'k');

hold on
semilogx(Q.sfprefdom,Q.DR.BWlog_2Cmplx,'g')
hold on
semilogx(Q.sfprefdom,Q.SI.BWlog_2Cmplx,'r')

%hold on
%loglog(sfpref_sim,BWVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)

xlabel('sfreq')
ylabel('log BW')
%xlim(([plotdom(1) plotdom(end)]))


%% Plot

plotdom = [1/32 1/16 1/8 1/4 1/2 1];

plotSizeVPrediction2(Q.meas.size,Q.SI_D.size,plotdom); %Measured size vs. scale invariant prediction
subplot(2,2,1)
hold on
hold on,loglog(Q.SI.size,Q.DR.size_2Cmplx,'g')

% id = find(~isnan(TCroAll.sfpref));
% plotSizeVPrediction2(Rsize(id),getSizefromSFpref(sf2sigModel,TCroAll.sfpref(id)),dom)

subplot(2,2,1)
hold on,
%loglog(Q.SI.size,Q.DR.size,'--r','LineWidth',2);
%loglog(.25./sfpref_sim,XsigVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)
set(gca,'XTick',[1/16 1/8 .25 .5 1 2]/2,'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
set(gca,'YTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
ylabel('RF width; deg (sigma)')
xlabel('predicted width; deg (sigma) 1/(4*sfreq)')

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
subplot(2,2,1)
hold on
loglog(Q.SI.BWlin,Q.DR.BWlin_2Cmplx,'g')

subplot(2,2,1)
hold on 
%loglog(2*Q.sfprefdom/pi,Q.DR.BWlin,'--r','LineWidth',2);
set(gca,'XTick',plotdom,'YTick',plotdom)
%loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 

subplot(2,2,2)
%set(gca,'XTickLabel',{'1/4'; '1/2'; '1'; '2'; '4'; '8'})
xlabel('bandwidth; c/deg (sigma)')
ylabel('# neurons')
xlim(log2([plotdom(1) plotdom(end)]))

subplot(2,2,3)
%set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi')
ylabel('# neurons')
xlim(log2([plotdom(1) plotdom(end)]))


% smearsig = sqrt(mean(log2(1+Q.meas.BWlin./Q.meas.sfpref).^2 - log2(1+2/pi).^2))
% smearsig_cycles = (2.^smearsig-1).*Q.sfprefdom
% smearsig = smearsig.*ones(size(Q.sfprefdom));
% bwE = (2.^sqrt( log2(1+2/pi).^2 + smearsig.^2 ) - 1) .* Q.sfprefdom;
% subplot(2,2,1)
% hold on 
% loglog(2*Q.sfprefdom/pi,bwE)


%Fit line in log BW domain, then convert back and put it in these axes 

BWlog = log2(Q.meas.BWlin./Q.meas.sfpref + 1);
[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),BWlog,[-.5 1]);

sigE = (2.^(xhat(1)*log2(Q.sfprefdom)+xhat(2)) - 1)  .*Q.sfprefdom;
subplot(2,2,1)
hold on
loglog(2*Q.sfprefdom/pi,sigE,'g')



% y = log2(Q.meas.BWlin./Q.meas.sfpref + 1);
% H = [log2(Q.meas.sfpref) ones(size(Q.meas.sfpref))];
%xhat = inv(H'*H)*H'*y;

% figure,plot(H(:,1),y,'.')
% hold on
% plot(Hhat(:,1),yhat)





[r p] = corrcoef(log2(Q.meas.sfpref),BWlog);
Hhat = [log2(Q.sfprefdom) ones(size(Q.sfprefdom))];
yhat =  Hhat*xhat(:);
figure,plot(log2(Q.meas.sfpref),BWlog,'.')
hold on
plot(Hhat(:,1),yhat)
ylabel('logBW')
xlabel('logSF')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])

%%
% 
% set(gca,'XTick',BWdom,'YTick',BWdom)
% loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
% xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 
% 
% 
% 
% plotSizeVPrediction2(Q.meas.sfpref,2.^Q.meas.BWlog,plotdom)
% subplot(2,2,1)
% hold on
% loglog(Q.SI.BWlog,Q.DR.BWlog_2Cmplx,'g')
% xlabel('sfpref')


%%
% figure
% subplot(2,2,1)
% hold on 
% loglog(2*sfprefdom_DLR/pi,BW_DLR/2,'--r','LineWidth',2);
% set(gca,'XTick',BWdom,'YTick',BWdom)
% loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
% xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 
% 
% subplot(2,2,2)
% %set(gca,'XTickLabel',{'1/4'; '1/2'; '1'; '2'; '4'; '8'})
% xlabel('bandwidth; c/deg (sigma)')
% ylabel('# neurons')
% xlim(log2([plotdom(1) plotdom(end)]))
% 
% subplot(2,2,3)
% %set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
% xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi')
% ylabel('# neurons')
% xlim(log2([plotdom(1) plotdom(end)]))





%%






%Get complex model from scale invariant pool
smearsigx = size_meas.^2 - size_SI
size_SIpool = sqrt(size_SI.^2 + smearsig.^2);  %1sig
BWlin_SIpool = 1./(size_SI*2*pi); %1sig
BWlog_SIpool = (log2(1-2/pi) - log2(1+2/pi)) * ones(size(sfprefdom))/2; %1sig
oriBW_SIpool = atan(BWlin_SI./sfprefdom)*180/pi;





%fit from Ringach data
sfprefdom_DLR = logspace(log10(1/16),log10(12),30);
sigE_DLR = getSizefromSFpref(sf2sigModel,sfprefdom_DLR);
BW_DLR = 2./(sigE_DLR*2*pi); %2sig bandwidth

Rsize = TCopAll.profileSize/2; %actual RF width (sigma)
Rsizepred = 1./(4*TCroAll.sfpref); %predicted width (sigma) from spatial frequency
dom = [1/32 1/16 1/8 1/4 1/2 1];

%%%%%Size%%%%%%%%%%%

id = find(~isnan(TCroAll.sfpref.*Rsize));
sfdom = logspace(log10(.25),log10(10),30);
[CmplxModel_RFsize smearsig] = getSmearModel(TCroAll.sfpref(id),Rsize(id),sfdom,sf2sigModel,'size');
smearsig
plotSizeVPrediction2(Rsize,Rsizepred,dom)
subplot(2,2,1)
hold on
hold on,loglog(1./(4*sfdom),CmplxModel_RFsize,'g')

% id = find(~isnan(TCroAll.sfpref));
% plotSizeVPrediction2(Rsize(id),getSizefromSFpref(sf2sigModel,TCroAll.sfpref(id)),dom)

subplot(2,2,1)
hold on,
loglog(.25./sfprefdom_DLR,sigE_DLR,'--r','LineWidth',2);
loglog(.25./sfpref_sim,XsigVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)
set(gca,'XTick',[1/16 1/8 .25 .5 1 2]/2,'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
set(gca,'YTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
ylabel('RF width; deg (sigma)')
xlabel('predicted width; deg (sigma) 1/(4*sfreq)')

subplot(2,2,2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('RF width; deg (sigma)')
ylabel('# neurons')
xlim(log2([dom(1) dom(end)]))

subplot(2,2,3)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted width; deg (sigma) 1/(4*sfreq)')
ylabel('RF width; deg (sigma)')
xlim(log2([dom(1) dom(end)]))

%%%%%Linear Bandwidth%%%%%%%%%%%
Rsize = TCroAll.sfBWLin/2; %actual BW (sigma)
Rsizepred = 2*TCroAll.sfpref/pi;%predicted bw (sig) from sf preference
%dom = logspace(log10(.25),log10(8),6);
dom = [1/4 1/2 1 2 4 8];

id = find(~isnan(TCroAll.sfpref.*Rsize));
sfdom = logspace(log10(.25),log10(10),30);
[CmplxModel_LinBW smearsig] = getSmearModel(TCroAll.sfpref(id),Rsize(id),sfdom,sf2sigModel,'BW');
smearsig
plotSizeVPrediction2(Rsize,Rsizepred,dom)
subplot(2,2,1)
hold on
hold on,loglog(2*sfdom/pi,CmplxModel_LinBW,'g')

subplot(2,2,1)
hold on 
loglog(2*sfprefdom_DLR/pi,BW_DLR/2,'--r','LineWidth',2);
set(gca,'XTick',BWdom,'YTick',BWdom)
loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 

subplot(2,2,2)
%set(gca,'XTickLabel',{'1/4'; '1/2'; '1'; '2'; '4'; '8'})
xlabel('bandwidth; c/deg (sigma)')
ylabel('# neurons')
xlim(log2([dom(1) dom(end)]))

subplot(2,2,3)
%set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi')
ylabel('# neurons')
xlim(log2([dom(1) dom(end)]))


%%%%%Log Bandwidth%%%%%%%%%%%
BWdom = [.5 1 2 4 8 16]/2;

%fit from Ringach data
sfprefdom_DLR = logspace(log10(1/4),log10(8),30);
sigE_DLR = getSizefromSFpref(sf2sigModel,sfprefdom_DLR);
BW_DLR = 2./(sigE_DLR*2*pi); %2sig bandwidth
BW_DLR = log2((sfprefdom_DLR+BW_DLR/2)./(sfprefdom_DLR)) + 2;

Rsize = TCroAll.sfBW; %actual BW (sigma)
Rsizepred = (1+2/pi)/(1-2/pi) * ones(1,length(Rsize));%predicted bw (sig) from sf preference
%dom = logspace(log10(.25),log10(8),6);
dom = log2([1/4 1/2 1 2 4 8]);

id = find(~isnan(TCroAll.sfpref.*Rsize));
sfdom = logspace(log10(.25),log10(10),30);
%[CmplxModel smearsig] = getSmearModel(TCroAll.sfpref(id),Rsize(id),sfdom,sf2sigModel,'BW');
smearsig
%plotSizeVPrediction2(Rsize,Rsizepred',dom)
figure,
subplot(2,2,1)
scatter(log2(TCroAll.sfpref),Rsize,'.k')
%hold on
%loglog(2*sfdom/pi,CmplxModel)

subplot(2,2,1)
hold on 
plot(log2(sfprefdom_DLR),BW_DLR,'--r','LineWidth',2);
set(gca,'XTick',BWdom,'YTick',BWdom)
plot(log2(sfpref_sim),BWVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 

subplot(2,2,2)
%set(gca,'XTickLabel',{'1/4'; '1/2'; '1'; '2'; '4'; '8'})
xlabel('bandwidth; c/deg (sigma)')
ylabel('# neurons')
xlim(log2([dom(1) dom(end)]))

subplot(2,2,3)
%set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi')
ylabel('# neurons')
xlim(log2([dom(1) dom(end)]))

%%
figure,
subplot(2,3,1)
semilogx(TCroAll.sfpref,TCroAll.F1F0,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency preference')
ylabel('F1/F0')
id = find(~isnan(TCroAll.sfpref.*TCroAll.F1F0));
[r p] = corrcoef(TCroAll.sfpref(id),TCroAll.F1F0(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

[F1F0model] = getF1F0vsSFmodel(Q.sfprefdom(:),xsigDom);
hold on, semilogx(Q.sfprefdom,F1F0model,'g')
xlim([.25 8])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})


subplot(2,3,2)
semilogx(TCroAll.sfBWLin/2,TCroAll.F1F0,'.k')
%hold on, semilogx(BWLinVsf_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency BW (sigma)')
ylabel('F1/F0')
id = find(~isnan(TCroAll.F1F0.*TCroAll.sfBWLin));
[r p] = corrcoef(TCroAll.F1F0(id),TCroAll.sfBWLin(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.SI.BWlin,F1F0model,'g')
xlim([.6 3])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})

[TCopAll.D TCopAll.Dnorm TCopAll.rOverlap] = OnOffOverlap(TC_op);

subplot(2,3,3)
plot(TCopAll.Dnorm,TCroAll.F1F0,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
%set(gca,'XTick',[.5 1 2 4 8])
xlabel('F1F0')
ylabel('Dnorm')
id = find(~isnan(TCroAll.F1F0.*TCopAll.Dnorm'));
[r p] = corrcoef(TCroAll.F1F0(id),TCopAll.Dnorm(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

subplot(2,3,4)
semilogx(TCroAll.sfpref,TCopAll.Dnorm,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency preference')
ylabel('Dnorm')
id = find(~isnan(log2(TCroAll.sfpref).*TCopAll.Dnorm'));
[r p] = corrcoef(log2(TCroAll.sfpref(id)),TCopAll.Dnorm(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

subplot(2,3,5)
semilogx(TCroAll.sfBWLin,TCopAll.Dnorm,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency BW (2sig)')
ylabel('Dnorm')
id = find(~isnan(log2(TCroAll.sfBWLin).*TCopAll.Dnorm'));
[r p] = corrcoef(log2(TCroAll.sfBWLin(id)),TCopAll.Dnorm(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

subplot(2,3,6)
plot(1./(2*TCroAll.sfpref),TCopAll.D,'.k')
hold on
plot([.1 1],[.1 1])
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
%set(gca,'XTick',[.5 1 2 4 8])
xlabel('1/(2sf)')
ylabel('D ')
id = find(~isnan(log2(TCroAll.sfpref).*TCopAll.D'));
[r p] = corrcoef(log2(TCroAll.sfpref(id)),TCopAll.D(id));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

DnormQuartiles = [prctile(TCopAll.Dnorm,0) prctile(TCopAll.Dnorm,25) prctile(TCopAll.Dnorm,50) prctile(TCopAll.Dnorm,75) prctile(TCopAll.Dnorm,100)]
F1F0Quartiles = [prctile(TCroAll.F1F0,0) prctile(TCroAll.F1F0,25) prctile(TCroAll.F1F0,50) prctile(TCroAll.F1F0,75) prctile(TCroAll.F1F0,100)]
%%
figure,plot(1./(2*TCroAll.sfpref),D,'.k')
hold on, plot([.01 1],[.01 1])

% subplot(1,3,3)
% semilogx(1./TCroAll.orisig,TCroAll.F1F0,'.k')
% hold on, semilogx(BWLinVsf_sim,orisigVsf,'-ob')
% %set(gca,'XTick',[.5 1 2 4 8])
% xlabel('orientation BW')
% ylabel('F1/F0')

%% Ori bandwidth vs. scale parameters

% aspectRatio = mean(Allparam(:,3)./Allparam(:,4));
% aspectRatio = 1;
% sfBWpred = 2*TCroAll.sfpref/pi; %predicted bw (sig) from sf preference
% sfBWpred_orth = sfBWpred/aspectRatio;
% oriBWpred = atan(sfBWpred_orth./TCroAll.sfpref)*180/pi;  %sig

%Just do the empirical measurement for the 5 bands
if strcmp(sf2sigModel,'NL sat')
    oriBWpred = [46 21 21 21 21]; %NL sat
elseif strcmp(sf2sigModel,'exp decay')
    oriBWpred = [37 20 20 20 20]; %exp dec
elseif strcmp(sf2sigModel,'classic')
    oriBWpred = [30 30 30 30 30]; %classic (scale invariant)
end

[sfprefdomI oriBWpredI] = getOriSigfromModel(sf2sigModel);


[sfprefdomI oriBW_scaleInv_pred] = getOriSigfromModel('classic');

%%
figure
subplot(2,1,1)
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));
loglog(TCroAll.sfpref,TCroAll.orisig,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),log2(TCroAll.orisig(id)));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 8.5])
ylim([8 45])

[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),log2(Q.meas.oriBW),[-.5 1]);
hold on
loglog(2.^dom,2.^ffit)


AR = 3;
%sfprefdom = logspace(log10(1/2),log10(8),30);
sizesig = getSizefromSFpref(sf2sigModel,sfprefdom_DLR);
BW = 1./(sizesig*2*pi); %1sig bandwidth
oripred = atan(BW./sfprefdom_DLR/AR)*180/pi;
hold on,
% plot(sfprefdom,oripred,'g')
% oripred = atan(1./(2*pi*.38*AR))*180/pi*ones(size(sfprefdom));

sfprefdom_DLR = logspace(log10(1/16),log10(12),30);
sigE_DLR = getSizefromSFpref(sf2sigModel,sfprefdom_DLR);
BW_DLR = 2./(sigE_DLR*2*pi); %2sig bandwidth

hold on,
plot(sfprefdom_DLR,oripred)

oripred_DLR = atan(BW_DLR./sfprefdom_DLR/AR/2)*180/pi;
hold on,
plot(sfprefdom_DLR,oripred_DLR)


%%
figure
subplot(2,1,1)
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));
loglog(TCroAll.sfpref,TCroAll.orisig,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),log2(TCroAll.orisig(id)));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 8.5])
ylim([8 45])


%hold on
%loglog(sfpref_sim,orisigVsf,'o-b','LineWidth',2,'MarkerSize',5)
hold on
loglog(sfprefdomI,oriBWpredI,'--r','LineWidth',2)
hold on
%loglog(sfprefdomI,oriBW_scaleInv_pred,'k')
loglog(sfprefdomI,oriBW_scaleInv_pred(1)*ones(1,length(sfprefdomI)),'k')


H = [log2(TCroAll.sfpref(id)) ones(size(TCroAll.sfpref(id)))];
pa = inv(H'*H)*H'*log2(TCroAll.orisig(id));

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'deg/octave'])
%axis square

%%
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));

sizepred = 1/(4*TCroAll.sfpref);
BWpred = 2*TCroAll.sfpref/pi; %2sig bandwidth

sfdum = [0.5 1 2 4 8];
BWpred = 2*TCroAll.sfpref/pi;
BWpred.^2 - TCroAll.sfpref

oriBWpred = atan(2/pi)*180/pi


oriBWpred = atan(1./(sfdum*signew))


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
%%

subplot(2,1,2)

BWsig = TCroAll.sfBWLin/2;
BWsig_sim = BWLinVsf_sim/2;
id = find(~isnan(BWsig.*TCroAll.orisig));
loglog(BWsig,TCroAll.orisig,'.k'), xlabel('SF bandwidth (sigma, c/deg)'), ylabel('ori bandwidth (sigma, deg)')
[r p] = corrcoef(log2(BWsig(id)),log2(TCroAll.orisig(id)));
set(gca,'XTick',[logspace(log10(.5),log10(4),4)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.5 4])
ylim([5 45])


%hold on
%loglog(BWsig_sim,orisigVsf,'o-b','LineWidth',2,'MarkerSize',5)
hold on
loglog(BWsig_sim,oriBWpredI,'--r','LineWidth',2)
hold on
loglog(BWsig_sim,oriBW_scaleInv_pred,'k')

H = [log2(BWsig(id)) ones(size(BWsig(id)))];
pa = inv(H'*H)*H'*TCroAll.orisig(id);

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'deg/octave'])
axis square


%% Compute pooling window in mm of cortical space

id = find(~isnan(TCroAll.sfpref.*TCopAll.profileSize));
xsig = getSizefromSFpref(sf2sigModel,TCroAll.sfpref(id));

1/(4*TCroAll.sfpref(id));
Rsizepred = xsig*2;

Rsize = TCopAll.profileSize(id); %actual RF width
%Rsize = TCopAll.ysize;
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
%id = 1:length(Rsize);
smearSig = sqrt(mean((Rsize(id)/2).^2 - (Rsizepred(id)/2).^2))
hist(real(smearSig))
xlabel('smear = sqrt(sigMeas^2 - sigPred^2)')
muSmear = 2*mean(smearSig); %degrees 2sigma
title(['Smear (2sig)=' num2str(muSmear)])

magFac = 1; %mm/deg
corticalSmear = muSmear*magFac; %mm 2sigma


%% Compute pooling window in cyc/deg of cortical space. i.e. the smearing function that creates wider sf tuning

id = find(~isnan(TCroAll.sfpref.*TCroAll.sfBWLin));
sigE = getSizefromSFpref(sf2sigModel,TCroAll.sfpref(id));

BWpred = 2./(sigE*2*pi); %2sigm Using Ringach data

BWpred = 2*2*TCroAll.sfpref(id)/(pi); %2sigm  Using standard model
BW = TCroAll.sfBWLin(id); %actual BW

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
id = 1:length(BW);
smearSig = sqrt(mean((BW(id)/2).^2 - (BWpred(id)/2).^2));
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
Dmax = 525;
Dbins = 0:75:Dmax;

e = 2;
SFOriPhaseClustering(PW_ro{2},Dbins)

sizeClustering(PW_op{2},Dbins)

sizeClustering(PWopAll,Dbins)
SFOriPhaseClustering(PWroAll,Dbins)

RetvPhaseClustering(PWroAll,PWopAll)
%%
figure, 
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));
loglog(TCroAll.sfpref,TCroAll.orisig,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),log2(TCroAll.orisig(id)));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
xlim([.2 8.5])

H = [log2(TCroAll.sfpref(id)) ones(size(TCroAll.sfpref(id)))];
pa = inv(H'*H)*H'*TCroAll.orisig(id)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'deg/octave'])

hold on

plot(sfpref,orisigVsf)

figure, 
AR = (TCopAll.xsize./TCopAll.ysize);
id = find(~isnan(AR.*TCroAll.sfpref));
loglog(TCroAll.sfpref,AR,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('Aspect ratio of RF')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),log(AR(id)));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
xlim([.2 8.5])

H = [log2(TCroAll.sfpref(id)) ones(size(AR(id)))];
pa = inv(H'*H)*H'*AR(id)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'octave/octave'])


%%
id = find(PWroAll.Dist>0 & PWroAll.Dist<50);
dori = abs(PWroAll.dori(id));
musf = sqrt(PWroAll.dsfpair(id,1).*PWroAll.dsfpair(id,2));
%musf = .5*(PWroAll.dsfpair(id,1)+PWroAll.dsfpair(id,2));

% an = 3;
% id = find(PW_ro{an}.Dist{1}>0 & PW_ro{an}.Dist{1}<100);
% dori = abs(PW_ro{an}.dori{1}(id));
% musf = sqrt(PW_ro{an}.dsfpair{1}(id,1).*PW_ro{an}.dsfpair{1}(id,2));
% %musf = .5*(PWroAll.dsfpair(id,1)+PWroAll.dsfpair(id,2));

figure, scatter(musf,dori,'.')

id = find(~isnan(dori(:).*musf(:)));
dori = dori(id);
musf = musf(id);
[r p] = corrcoef(musf,dori)


%% Compare size to predicted size from spatial frequency preference

sizeDom = logspace(log10(1/(2*getparam('max_sf'))),log10(1/(2*getparam('min_sf'))),2*getparam('n_sfreq'))/2;
BWdom = [.5 1 2 4 8 16]/2;

%fit from Ringach data
sfprefdom_DLR = logspace(log10(1/16),log10(12),30);


sigE_DLR = getSizefromSFpref(sf2sigModel,sfprefdom_DLR);

BW_DLR = 2./(sigE_DLR*2*pi); %2sig bandwidth

figure, 

%%%%RF Width vs. prediction of width, using spatial frequency preference
Rsize = TCopAll.profileSize/2; %actual RF width (sigma)
%Rsize = TCopAll.ysize/2;
Rsizepred = 1./(4*TCroAll.sfpref); %predicted width (sigma) from spatial frequency
subplot(3,2,1)
plotSizeVPrediction_scatter(Rsize,Rsizepred,sizeDom)
hold on,
loglog(.25./sfprefdom_DLR,sigE_DLR,'--r','LineWidth',2);

loglog(.25./sfpref_sim,XsigVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)
set(gca,'XTick',[1/16 1/8 .25 .5 1 2]/2,'YTick',[1/16 1/8 .25 .5 1 2]/2)
set(gca,'XTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
set(gca,'YTickLabel',{'1/32'; '1/16'; '1/8'; '1/4'; '1/2'; '1'})
ylabel('RF width; deg (sigma)')
xlabel('predicted width; deg (sigma) 1/(4*sfreq)')
%xlim([.05 1.5]), ylim([.05 1.5])

subplot(3,2,2), 
plotSizeVPrediction_hist(Rsize,Rsizepred)

%%%%%%SF bandwidth vs. prediction of bandwidth, using spatial frequency preference
subplot(3,2,3)
Rsize = TCroAll.sfBWLin/2; %actual BW (sigma)
Rsizepred = 2*TCroAll.sfpref/pi;%predicted bw (sig) from sf preference
plotSizeVPrediction_scatter(Rsize,Rsizepred,BWdom)
xlabel('predicted Bandwidth; c/deg (sigma): 2sfreq/pi'), ylabel('bandwidth; c/deg (sigma)') 
hold on
loglog(2*sfprefdom_DLR/pi,BW_DLR/2,'--r','LineWidth',2);
set(gca,'XTick',BWdom,'YTick',BWdom)

loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
subplot(3,2,4), 
plotSizeVPrediction_hist(Rsize,Rsizepred)

%%%%%%Ori bandwidth vs. spatial frequency preference
aspectRatio = 2;
sfBWpred = 2*TCroAll.sfpref/pi; %predicted bw (sig) from sf preference
oriBWpred = atan(2/pi/aspectRatio)*180/pi;  %sig


subplot(3,2,5)
id = find(~isnan(TCroAll.sfpref.*TCroAll.orisig));
loglog(TCroAll.sfpref,TCroAll.orisig,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(TCroAll.sfpref(id)),log2(TCroAll.orisig(id)));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 8.5])
ylim([5 45])


hold on
loglog(sfpref_sim,orisigVsf,'o-b','LineWidth',2,'MarkerSize',5)
hold on
plot([.25 7],[oriBWpred oriBWpred],'k')

H = [log2(TCroAll.sfpref(id)) ones(size(TCroAll.sfpref(id)))];
pa = inv(H'*H)*H'*TCroAll.orisig(id);

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(pa(1)) 'deg/octave'])
axis square

