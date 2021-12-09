function Ggetrevcorrkernel2(cellMat,trialdom,hh)

%Ian Nauhaus

%Keep cellMat as an input because it filters it, and I don't want to
%create a new variable of the same size

%2 computes the standard deviation as well. 

%4 builds the kernel baased on an ARMA model

global ACQinfo Analyzer cellS maskS G_RChandles


DC = 0;
ARflag = 1;
polyorder = 1;
blankNorm = get(G_RChandles.blankNorm,'value');

%%%Get rid of the glia:

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);
Ncell = length(nID);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)
[dum id0] = min(abs(taudom-0));  
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

if ~exist('frate')
    frate = 60.0044;
end

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file)

[domains seqs] = getSeqInfo(trialdom);


%%%%%%%%%%%%%%%%%%%

oridom = domains{trialdom(1)}.oridom;
sfdom = domains{trialdom(1)}.sfdom;
phasedom = domains{trialdom(1)}.phasedom;
colordom = domains{trialdom(1)}.colordom;

Ncell = length(nID);

cellS.muTime = cell(1,Ncell);
cellS.sigTime = cell(1,Ncell);
cellS.muBase = zeros(1,Ncell);
cellS.sigBase = zeros(1,Ncell);
cellS.muTimeBlank = cell(1,Ncell);
cellS.sigTimeBlank = cell(1,Ncell);

%paramP = length(oridom)*length(sfdom)*length(phasedom);
paramP = length(oridom)*length(sfdom);

if ~ARflag
    tauP = round(500/acqPeriod)+1;
else
    tauP = round(260/acqPeriod)+1;
end
%tauP = 40;
taudom = 0:acqPeriod:((tauP-1)*acqPeriod);
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end


%Preallocate
Precountmat = getGratingPresentationCounts(trialdom);
for ori = 1:length(oridom)
    for sf = 1:length(sfdom)
        for phase = 1:length(phasedom)
            for color = 1:length(colordom)
                kernPre{ori,sf,phase,color} = zeros(Precountmat(ori,sf,phase,color),length(taudom));
            end
        end
    end
end

NT = getnotrials;
tcoursedum = 0;
figure
for p = 1:Ncell
tic
    pID = nID(p);

    [idcelly idcellx] = find(masklabel == celldom(pID));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom)); 
    countmatblank{p} = 0;
    
    %kern{p} = NaN*zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));    
    %kernblank{p} = zeros(1,length(taudom));


    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);

        hper = gethper(cond);
        Tupdate = Tf*hper;

        %         dsync = diff(cellS.synctimes{cond,rep});
        %         if any(abs(dsync-Tupdate/1000)>.100)
        %             'Warning: syncs may be messed up'
        %         end

       
         tcourse = squeeze(cellMat{cond}(pID,:,rep));
        
        id = find(isnan(tcourse));
        tcourse(id) = nanmean(tcourse);

        if ~isnan(sum(tcourse))  %some cells go out of the field of view

            y = processTcourse(tcourse,hh,1,acqPeriod);              
%             id = find(tcourse<prctile(tcourse,50) & tcourse>prctile(tcourse,0));
%             tcourse = tcourse-median(tcourse); 
%             tcourse = tcourse/std(tcourse(id));
            

            tdom = (0:length(y)-1)*acqPeriod;
            %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
            tdom_pix = tdom + tau_xy;
 
            %Deconvolution:
            %Ctrx = exp(-tdom/900);
            %tcourse = ifft(fft(tcourse)./fft(Ctrx));


            %idbase = find(tdom_pix<cellS.synctimes{cond,rep}(1)*1000 | tdom_pix>(cellS.synctimes{cond,rep}(end)*1000+500));
            %bLine = mean(tcourse(idbase));

            %tcourse = (tcourse-bLine)/bLine;

            
            %HHdum = zeros(length(y),length(oridom)*length(sfdom)*length(phasedom));
            HHdum = zeros(length(y),length(oridom)*length(sfdom));
            
            pid = 1;
            for ori = 1:length(oridom)
                for sf = 1:length(sfdom)
                    %for phase = 1:length(phasedom)
                        
                        %id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf )& seqs{T}.phaseseq == phasedom(phase));
                        id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                        %id = find(seqs{T}.oriseq == oridom(ori));
                        
                        if ~isempty(id)
                            
                            stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)
                            
                            for i = 1:length(stimes)
                                [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                                HHdum(idx1,pid) = 1;
                            end
                        end
                        pid = pid+1;
                   % end
                end
            end
            
            %HHdum = zscore(HHdum);
            HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
            HHdum = [.0*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
            for z = tauP:length(HHdum(:,1))
                chunk = squeeze(HHdum(z-tauP+1:z,:))';
                chunk = chunk(:)';
                HH(z-tauP+1,:) = chunk;
            end
            
            if DC
                HH = [HH ones(length(HH(:,1)),1)];
            end
            if ARflag
                
                %dydt = [(y(2:end)-y(1:end-1)) 0];
                %dydt = [0 (y(3:end)-y(1:end-2))/2 0];
                %HH = [HH dydt'];
                
                HH = [HH [0; y(1:end-1)']];
            end
            
            covMat{p} = covMat{p} + HH'*HH;
            
            xCorr{p} = HH'*y(:) + xCorr{p};
            
%             id = find(isnan(seqs{T}.oriseq));
%             stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times
            %stimes = cellS.synctimes{cond,rep}(id)*1000;

            %stimes(find(stimes>stimesW(2) | stimes<stimesW(1))) = [];

% 
%             for i = 1:length(stimes)
% 
%                 [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
%                 if ~isempty(idx1)
% 
%                     idx1 = idx1(1);
%                     tpiece = idx1:idx1+length(taudom)-1;
%                     tbase = (idx1+id0-1)-1:(idx1+id0-1)-1;  %Use time zero as base
% 
%                     if tpiece(1)>0 & tpiece(end)<length(tcourse)
%                         base = mean(squeeze(tcourse(tbase)));
%                         tcoursePiece = squeeze(tcourse(tpiece))';
% 
%                         %tcoursePiece = (tcoursePiece-base);
%                         kernblank{p}(countmatblank{p}+1,:) = tcoursePiece;
%                         countmatblank{p} = countmatblank{p} + 1;
%                     end
%                 end
%             end

        end
    end
    
    
    
        
    if ~ARflag
        covMat{p} = covMat{p}.*eye(size(covMat{p}));
    end
    
    params{p} = (covMat{p})\xCorr{p};    
    
    %params{p} = xCorr{p};
    
    paramsdum = params{p};
    if ARflag
        alpha = params{p}(end);
        tauConst(p) = acqPeriod/-log(alpha);
        paramsdum = params{p}(1:end-1);
        xCorr{p} = xCorr{p}(1:end-1);
    end
    if DC
        xCorr{p} = xCorr{p}(1:end-1);
        paramsdum = paramsdum(1:end-1);  %Get rid of DC shift
    end

    kerndum = reshape(paramsdum,paramP,tauP);

    for i = 1:length(kerndum(1,:)) %loop each time point
        dum = kerndum(:,i);
        %kern{p}(:,:,:,i) = reshape(dum,[length(sfdom) length(phasedom) length(oridom)]);  %first dimension will be ori
        kern{p}(:,:,:,i) = reshape(dum,[length(sfdom) length(oridom)]);  %first dimension will be ori
    end
    
 
    
    
% 
%     id = find(countmat{p} == 0);
%     if ~isempty(id)
%         'Some stimuli never shown.'
%     end
% 
%     if sum(countmatblank{1})
%         cellS.muTimeBlank{p} = mean(kernblank{p});
%         cellS.sigTimeBlank{p} = std(kernblank{p})/sqrt(length(kernblank{p})); %standard error
%     end
% 
%     cellS.kernAll = kern;
% 
%     kern = [];
% 
%     for ori = 1:length(oridom)
%         for sf = 1:length(sfdom)
%             for phase = 1:length(phasedom)
%                 for color = 1:length(colordom)
%                     if ~isempty(cellS.kernAll{p}{ori,sf,phase,color}) %if it didn't show this grating
%                         N = length(cellS.kernAll{p}{ori,sf,phase,color}(:,1)); %Number of presentations
% 
%                         if blankNorm
%                             cellS.muTime{p}(ori,sf,phase,color,:) = squeeze(mean(cellS.kernAll{p}{ori,sf,phase,color},1) - cellS.muTimeBlank{p});% / mean(cellS.muTimeBlank{p});
%                             cellS.sigTime{p}(ori,sf,phase,color,:) = std(cellS.kernAll{p}{ori,sf,phase,color},[],1)/sqrt(N);% / mean(cellS.muTimeBlank{p});
%                         else
%                             cellS.muTime{p}(ori,sf,phase,color,:) = squeeze(mean(cellS.kernAll{p}{ori,sf,phase,color},1));
%                             cellS.sigTime{p}(ori,sf,phase,color,:) = std(cellS.kernAll{p}{ori,sf,phase,color},[],1)/sqrt(N);
%                         end
%                     end
%                 end
%             end
%         end
%    end

    %Subtract mean from each time point
    % for z = 1:length(cellS.muTime{p}(1,1,1,1,:))
    %
    %     dum = cellS.muTime{p}(:,:,:,:,z);
    %     cellS.muTime{p}(:,:,:,:,z) = cellS.muTime{p}(:,:,:,:,z) - mean(dum(:));
    %
    % end

    %Base is used for data selection later...
    %     dum = cellS.muTime{p}(:,:,:,:,1);
    %     cellS.muBase(p) = mean(dum(:));  %average response at time zero
    %     dum = cellS.sigTime{p}(:,:,:,:,1);
    %     cellS.sigBase(p) = mean(dum(:));  %average sigma response at time zero


    %kernplot = squeeze(mean(mean(kern{p},1),2));
    kernplot = squeeze(mean(kern{p},1));
    
   
    
%     kernplot = cellS.muTime{p}(:,:,:,:,:);
%     kernplot = mean(kernplot,3);  %mean across phase
%     kernplot = mean(kernplot,4); %mean across color

    %kernblank = kern{p}(end,:);

    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    
    %blankmat = ones(length(kernplot(:,1)),1)*kernblank{p}(:)';

    %kernplot = (kernplot-blankmat);
    %kernplot = [kernplot; kernblank{p}(:)'];

    %imagesc(round(taudom),sfdom,real(fliplr(kernplot)))

    %dum = dum./(ones(size(dum,1),1)*max(dum));
    imagesc((kernplot))
    
    for q = 1:size(kernplot,2)
        [param ffit_ori varacc_ori(p,q) sigma] = Gaussfit(oridom,kernplot(:,q)',1);
        Gsig(p,q) = param(2);
    end
    
    %plot(cellS.muTimeBlank{p})

    drawnow
    toc
end

id = find(varacc_ori<.8);
Gsigdum = Gsig;
Gsigdum(id) = NaN;

figure,scatter(Gsigdum(:,3),Gsigdum(:,2)), hold on, plot([0 100],[0 100],'r'),xlim([0 100]),ylim([0 100])

'hi'



