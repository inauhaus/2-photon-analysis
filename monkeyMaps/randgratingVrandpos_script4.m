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

global alpha
alpha = 4;

sf2sigModel = 'NL sat';
%sf2sigModel = 'exp decay';
sf2sigModel = 'classic'

sfBands = linspace(1,5,5);
if strcmp(sf2sigModel,'NL sat') 
    xsigDom = .2;
    sfsigDom = 1;
elseif strcmp(sf2sigModel,'classic') 
    xsigDom = .2;
    sfsigDom = .75;
end
%sfsigDom = .5
%[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation_map(sfBands,xsigDom,sfsigDom,sf2sigModel);

[XsigVsf_sim BWLinVsf_sim BWVsf_sim sfpref_sim F1vsBand_sim] = V1_ComplexCell_generation(sfBands,xsigDom,sfsigDom,sf2sigModel);
%[XsigVsf_sim BWLinVsf_sim sfpref_sim F1vsBand_sim orisigVsf] = V1_ComplexCell_generation_2Dx(sfBands,xsigDom,sfsigDom,sf2sigModel);

BWLinVsf_sim = BWLinVsf_sim*2; %Make it 2sigma to be consistant with stuff below.
%sfpref_sim = sfBands


%% Set up layers of the model, starting with the SF domain 

global alpha
alpha = pi;


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
Q.SI.size = 1./(alpha*Q.sfprefdom);  %1sig
Q.SI.BWlin = 1./(Q.SI.size*2*pi); %1sig
Q.SI.BWlog = log2(1+2/pi) * ones(size(Q.sfprefdom)); %1sig
Q.SI.oriBW = atan(Q.SI.BWlin./Q.sfprefdom)*180/pi;
%domain of measured SF preference
Q.SI_D.size = 1./(alpha*Q.meas.sfpref);  %1sig
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


Q.SI_D.size_2Cmplx = sqrt(Q.SI_D.size.^2 + Q.SIsmearXsig.^2); %Model fit of size, using scale invariance as inputs
Q.DR_D.size_2Cmplx = sqrt(Q.DR_D.size.^2 + Q.DRsmearXsig.^2); %Model fit of size, using simple cells as inputs

Q.SI_D.BWlin_2Cmplx = sqrt(Q.SI_D.BWlin.^2 + Q.SIsmearSFsig.^2);  %Model fit of linear SF bw, using scale invariance as inputs
Q.DR_D.BWlin_2Cmplx = sqrt(Q.DR_D.BWlin.^2 + Q.DRsmearSFsig.^2);  %Model fit of linear SF bw, using simple cells as inputs

Q.SI_D.BWlog_2Cmplx = log2((Q.SI_D.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref));  %Model fit of log SF bw, using scale invariance as inputs
Q.DR_D.BWlog_2Cmplx = log2((Q.DR_D.BWlin_2Cmplx+Q.meas.sfpref) ./ (Q.meas.sfpref)); %Model fit of log SF bw, using simple cells as input

%%
%%%%%Log Bandwidth%%%%%%%%%%%
% plotdom = [1/2 1 2 4];
% figure
% subplot(2,2,1)
% semilogx(Q.meas.sfpref,Q.meas.BWlog,'.k');
% hold on
% semilogx(Q.sfprefdom,Q.DR.BWlog,'--k');
% hold on
% semilogx(Q.sfprefdom,Q.SI.BWlog,'k');
% 
% hold on
% semilogx(Q.sfprefdom,Q.DR.BWlog_2Cmplx,'g')
% hold on
% semilogx(Q.sfprefdom,Q.SI.BWlog_2Cmplx,'r')
% 
% %hold on
% %loglog(sfpref_sim,BWVsf_sim,'o-b','LineWidth',2,'MarkerSize',5)
% 
% xlabel('sfreq')
% ylabel('log BW')
% %xlim(([plotdom(1) plotdom(end)]))


%% Plot

plotdom = [1/32 1/16 1/8 1/4 1/2 1];

plotSizeVPrediction2(Q.meas.size,Q.SI_D.size,plotdom); %Measured size vs. scale invariant prediction
subplot(2,2,1)
hold on
hold on,loglog(Q.DR.size,Q.DR.size_2Cmplx,'g')

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

[sig_x sig_ori] = SFxyPoolingModel(Q.sfprefdom);

%%%%%Linear Bandwidth%%%%%%%%%%%
plotdom = [1/4 1/2 1 2 4 8];

plotSizeVPrediction2(Q.meas.BWlin,Q.SI_D.BWlin,plotdom)
subplot(2,2,1)
hold on
loglog(Q.SI.BWlin,Q.SI.BWlin_2Cmplx,'g') 
hold on
loglog(Q.SI.BWlin,sig_x(1,:),'k') 


subplot(2,2,1)
hold on 
%loglog(2*Q.sfprefdom/pi,Q.DR.BWlin,'--r','LineWidth',2);
set(gca,'XTick',plotdom,'YTick',plotdom)
%loglog(2*sfpref_sim/pi,BWLinVsf_sim/2,'o-b','LineWidth',2,'MarkerSize',5)
xlabel('predicted Bandwidth; c/deg (sigma): sfreq/(alpha*pi*2)'), ylabel('bandwidth; c/deg (sigma)') 

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


%Fit line in log BW domain, then convert back and put it in these axes 

BWlog = log2(Q.meas.BWlin./Q.meas.sfpref + 1); %Equation A
[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),BWlog,[-1 1]);

sigE = (2.^(xhat(1)*log2(Q.sfprefdom)+xhat(2)) - 1)  .*Q.sfprefdom; %algebraic manipulation of Eq. A above gives BW(sfpref)
subplot(2,2,1)
hold on
plot(Q.DR.BWlin,sigE,'k')



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
idBP = find(TCroAll.LPness<=.61 & ~isnan(TCroAll.sfpref.*TCroAll.F1F0.*TCroAll.sfBWLin) & ~isnan(TCopAll.Dnorm(:))); %get bandpass population

DnormQuartiles = [prctile(TCopAll.Dnorm,0) prctile(TCopAll.Dnorm,25) prctile(TCopAll.Dnorm,50) prctile(TCopAll.Dnorm,75) prctile(TCopAll.Dnorm,100)]
F1F0Quartiles = [prctile(TCroAll.F1F0,0) prctile(TCroAll.F1F0,25) prctile(TCroAll.F1F0,50) prctile(TCroAll.F1F0,75) prctile(TCroAll.F1F0,100)]

%DnormQuartiles = [prctile(TCopAll.Dnorm(idBP),0) prctile(TCopAll.Dnorm(idBP),25) prctile(TCopAll.Dnorm(idBP),50) prctile(TCopAll.Dnorm(idBP),75) prctile(TCopAll.Dnorm(idBP),100)]
F1F0Quartiles = [prctile(TCroAll.F1F0(idBP),0) prctile(TCroAll.F1F0(idBP),25) prctile(TCroAll.F1F0(idBP),50) prctile(TCroAll.F1F0(idBP),75) prctile(TCroAll.F1F0(idBP),100)]

[r p] = corrcoef(TCroAll.F1F0(idBP),TCopAll.Dnorm(idBP))

figure,
subplot(2,3,1)
semilogx(TCroAll.sfpref(idBP),TCroAll.F1F0(idBP),'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
xlim([.25 8])
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency preference')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.sfpref(idBP),TCroAll.F1F0(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

[F1F0model] = getF1F0vsSFmodel(Q.sfprefdom(:),xsigDom);
hold on, semilogx(Q.sfprefdom,F1F0model,'g')
xlim([.25 8])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})

subplot(2,3,2)
semilogx(TCroAll.sfBWLin(idBP)/2,TCroAll.F1F0(idBP),'.k')
xlim([.5 4])
%hold on, semilogx(BWLinVsf_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency BW (sigma)')
ylabel('F1/F0')
[r p] = corrcoef(TCroAll.F1F0(idBP),TCroAll.sfBWLin(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

hold on, semilogx(Q.SI.BWlin_2Cmplx,F1F0model,'g')
%hold on, semilogx(Q.SI.BWlin,F1F0model,'g')
%xlim([.6 3])
set(gca,'YTick',[0 pi/2 pi],'YTickLabels',{'0','pi/2','pi'})

[TCopAll.D TCopAll.Dnorm TCopAll.rOverlap] = OnOffOverlap(TC_op);

subplot(2,3,3)
plot(TCopAll.Dnorm(idBP),TCroAll.F1F0(idBP),'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
%set(gca,'XTick',[.5 1 2 4 8])
xlabel('F1F0')
ylabel('Dnorm')
%id = find(~isnan(TCroAll.F1F0.*TCopAll.Dnorm'));
[r p] = corrcoef(TCroAll.F1F0(idBP),TCopAll.Dnorm(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

subplot(2,3,4)
semilogx(TCroAll.sfpref,TCopAll.Dnorm,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency preference')
ylabel('Dnorm')
id = find(~isnan(log2(TCroAll.sfpref).*TCopAll.Dnorm'));
[r p] = corrcoef(log2(TCroAll.sfpref(idBP)),TCopAll.Dnorm(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])

subplot(2,3,5)
semilogx(TCroAll.sfBWLin,TCopAll.Dnorm,'.k')
%hold on, semilogx(sfpref_sim,F1vsBand_sim,'-ob')
set(gca,'XTick',[.5 1 2 4 8])
xlabel('spatial frequency BW (2sig)')
ylabel('Dnorm')
id = find(~isnan(log2(TCroAll.sfBWLin).*TCopAll.Dnorm'));
[r p] = corrcoef(log2(TCroAll.sfBWLin(idBP)),TCopAll.Dnorm(idBP));
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
[r p] = corrcoef(log2(TCroAll.sfpref(idBP)),TCopAll.D(idBP));
title(['r=' num2str(r(1,2)) '; p=' num2str(p(1,2))])


%%
figure,plot(1./(2*TCroAll.sfpref),TCopAll.D,'.k')
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
subplot(3,1,1)
semilogx(Q.meas.sfpref,Q.meas.oriBW,'.k'), xlabel('spatial frequency preference (c/deg)'), ylabel('orientation bandwidth (sigma, deg)')
[r p] = corrcoef(log2(Q.meas.sfpref),log2(Q.meas.oriBW));
set(gca,'XTick',[logspace(log10(.25),log10(8),6)])
set(gca,'YTick',round([logspace(log10(5),log10(45),6)]))
xlim([.2 8.5])
ylim([8 45])

[xhat ffit dom varaccount] = LinetotalLS(log2(Q.meas.sfpref),log2(Q.meas.oriBW),[-.5 1]);
hold on
semilogx((2.^dom),2.^ffit)

title(['r = ' num2str(r(1,2)) '; p = ' num2str(p(1,2)) '; slope = ' num2str(xhat(1)) ])


AR = 1.7;
oripred_SA = atan(Q.SI.BWlin./Q.sfprefdom/AR)*180/pi; 
oripred = atan(Q.SI.BWlin_2Cmplx./Q.sfprefdom/AR)*180/pi - 5;
oripred_D = atan(Q.SI_D.BWlin_2Cmplx./Q.meas.sfpref/AR)*180/pi;  %Same thing but using a different SF domain

hold on,
loglog(Q.sfprefdom,oripred,'r')
hold on,
loglog(Q.sfprefdom,oripred_SA,'k')

[sig_x sig_ori] = SFxyPoolingModel(Q.sfprefdom,0.75+.2);
hold on
loglog(Q.sfprefdom,sig_ori(1,:),'k') 
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

muori = 
cos()

%% Assessing response amplitude to bars vs. SF tuning



idLow = (find(TCroAll.fhi<2.5));
idMid = (find(TCroAll.fhi>2.5 & TCroAll.flo<2.5));
idHi = (find(TCroAll.flo>2.5));

Ncells = length(~isnan(TCroAll.fhi));

% LowMedian = nanmedian(TCopAll.amp(idLow))
% HighMedian = nanmedian(TCopAll.amp(idnLow))
% HighMedian/LowMedian
% [h p] = ttest2(TCopAll.amp(idLow),TCopAll.amp(idnLow))


for i = 1:length(TCopAll.rawProfile)
    amp(i) = max(TCopAll.rawProfile{i});
end

LowMedian = nanmedian(amp(idLow))
MidMedian = nanmedian(amp(idMid))
HighMedian = nanmean(amp(idHi))
PercAttenuation_Hi = HighMedian/MidMedian
PercAttenuation_Lo = LowMedian/MidMedian

[h p] = ttest2(amp(idLow),amp(idMid))
[h p] = ttest2(amp(idHi),amp(idMid))

id = find(~isnan(TCroAll.sfpref.*TCopAll.amp))
[r p] = corrcoef(TCroAll.sfpref(id),TCopAll.amp(id))





