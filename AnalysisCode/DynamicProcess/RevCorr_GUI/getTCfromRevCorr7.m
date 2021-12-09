function getTCfromRevCorr7(DSflag)

%Ian Nauhaus

%3 is from 2, and it combines the spatial frequency and orientation loops

%5 allows for downsampling of the kernel (e.g. with color experiments)

global cellS DM TC MK idExamp kernC kernSigC G_RChandles ACQinfo


trialdom = 1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];
countmat = getGratingPresentationCounts(trialdom);

%Downsample kernel and domains
if DSflag
    oriDS = 1;
    sfDS = 0;
    phaseDS = 0;
    for p = 1:length(cellS.muTime)
        
        kernsum = zeros(size(cellS.muTime{p}));
        kernsigsum = zeros(size(cellS.muTime{p}));
        for i = 1:size(kernsum,5)
            kernsum(:,:,:,:,i) = cellS.muTime{p}(:,:,:,:,i).*countmat;  %now it is the sum, instead of the mean
            kernsigsum(:,:,:,:,i) = cellS.sigTime{p}(:,:,:,:,i).*countmat;  %now it is the sum, instead of the mean
        end
        
        %Downsample orientation
        if oriDS
            kernsum = kernsum(1:2:end,:,:,:,:) + kernsum(2:2:end,:,:,:,:);
            kernsigsum = kernsigsum(1:2:end,:,:,:,:) + kernsigsum(2:2:end,:,:,:,:);
            cdum = countmat(1:2:end,:,:,:) + countmat(2:2:end,:,:,:);
        else
            cdum = countmat;
        end
        
        %Downsample spatial frequency
        if sfDS
            if ~rem(size(kernsum,2),2)
                kernsum = kernsum(:,1:2:end,:,:,:) + kernsum(:,2:2:end,:,:,:);
                kernsigsum = kernsigsum(:,1:2:end,:,:,:) + kernsigsum(:,2:2:end,:,:,:);
                cdum = cdum(:,1:2:end,:,:) + cdum(:,2:2:end,:,:);
            else
                kernsum = kernsum(:,1:2:end-1,:,:,:) + kernsum(:,2:2:end,:,:,:);  %throw away the last SF
                kernsigsum = kernsigsum(:,1:2:end-1,:,:,:) + kernsigsum(:,2:2:end,:,:,:);
                cdum = cdum(:,1:2:end-1,:,:) + cdum(:,2:2:end,:,:);
            end
        end

        %Downsample spatial phase
        if phaseDS
            kernsum = sum(kernsum,3);
            kernsigsum = sum(kernsigsum,3);
            cdum = sum(cdum,3);
        end
        
        
        %Renormalize
        muTime{p} = zeros(size(kernsum));
        sigTime{p} = zeros(size(kernsigsum));
        for i = 1:size(muTime{p},5);
            muTime{p}(:,:,:,:,i) = kernsum(:,:,:,:,i)./cdum;  %now it is the mean, instead of the sum
            sigTime{p}(:,:,:,:,i) = kernsigsum(:,:,:,:,i)./cdum;  %now it is the mean, instead of the sum
        end
        
        
    end
    
    if oriDS
        DM.oridom = angle(exp(1i*2*DM.oridom(1:2:end)*pi/180) + exp(1i*2*DM.oridom(2:2:end)*pi/180))/2*180/pi;
        id = find(DM.oridom<0);
        DM.oridom(id) = DM.oridom(id)+180;
        DM.oridom = DM.oridom-min(DM.oridom);  %This is cheating
    end
    
    if sfDS
        if ~rem(size(kernsum,2),2)
            DM.sfdom = (DM.sfdom(1:2:end) + DM.sfdom(2:2:end))/2;
        else
            DM.sfdom = (DM.sfdom(1:2:end-1) + DM.sfdom(2:2:end))/2;
        end
    end
    
    if phaseDS
        DM.phasedom = 0;
    end


    'Downsampling because not enough presentations'
    
else
    muTime = cellS.muTime;
    sigTime = cellS.sigTime;
end


TC = struct;

kernC = cell(1,length(DM.colordom));

for i = 1:length(DM.colordom)
    for p = 1:length(cellS.muTime)
        kernC{i}{p} = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernSigC{i}{p} = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        
        %Need to preserve dimension, even if it is only one element (i.e. one orientation)
        kernC{i}{p}(:,:,:,:) = squeeze(muTime{p}(:,:,:,i,:));
        kernSigC{i}{p}(:,:,:,:) = squeeze(sigTime{p}(:,:,:,i,:));
    end    
end


%t2 = getparam('h_per')*10+150; %Keep it restricted to the "onset response" ??
t2 = 400;

delayWin = [100 t2];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-DM.taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-DM.taudom));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions

tdom = DM.taudom(1:end);

kernsmooth = getSmoother(30,10,.2,tdom,DM.oridom,DM.sfdom); %for establising time-to-peak

%Smoother before taking ori curve
%orismoother = getSmoother(10,5,.1,DM.taudom,DM.oridom,DM.sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother(10,5,.01,tdom,DM.oridom,DM.sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(30,10,.05,tdom,DM.oridom,DM.sfdom);  %this is used to determine a "significant" response, 
                                                                %so keep smoothing minimal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get blank response
ht = zeros(1,length(DM.taudom));
ht(1:3) = [.3 1 .3]; ht = ht/sum(ht);
for c = 1:length(DM.colordom)

    for p = 1:MK.Ncell
        if c == 1
            muBlank{p} = 0;
            sigBlank{p} = 0;
        end
       if getparam('blankProb')==0
            kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
            kernplot(:,:,:,:) = kernC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
                        
            kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
            muBlankdum = squeeze(mean(kernplot(:,end,:),1)); %average over ori; use last spatial frequency for the blank
            
            kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
            kernplot(:,:,:,:) = kernSigC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
            kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
            sigBlankdum = squeeze(mean(kernplot(:,end,:),1)); %average over ori; use last spatial frequency for the blank
        else
            muBlankdum = cellS.muTimeBlank{p};
            sigBlankdum = cellS.sigTimeBlank{p};
        end

        muBlankdum = ifft(fft(muBlankdum(:)).*abs(fft(ht(:))));
        muBlank{p} = muBlank{p} + muBlankdum/length(DM.colordom); %average across color
        
        sigBlankdum = ifft(fft(sigBlankdum(:)).*abs(fft(ht(:))));
        sigBlank{p} = sigBlank{p} + sigBlankdum/length(DM.colordom); %average across color
    end
    
end

%subtract the blank
for c = 1:length(DM.colordom)
    for p = 1:MK.Ncell
        
        [dum idz] = min(abs(DM.taudom-0));
        
        base = kernC{c}{p}(:,:,:,idz);
        base = mean(base(:));
        
        muBlankdum = muBlank{p}(:) - muBlank{p}(idz) + base; %give them the same value at t = 0;
        for or = 1:length(DM.oridom)
            for sfr = 1:length(DM.sfdom)
                for phs = 1:length(DM.phasedom)
                    kernC{c}{p}(or,sfr,phs,:) = squeeze(kernC{c}{p}(or,sfr,phs,:)) - muBlankdum;
                end
            end
        end
    end
end


if length(DM.colordom) == 3
    colordomdum = 1:4;
else 
    colordomdum = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1flag = 0;
if f1flag
    phasekern = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
    for i = 1:length(DM.phasedom)        
        phasekern(:,:,i,:) = exp(1i*DM.phasedom(i)*pi/180);        
    end
end

colorid = {'r','g','b','k'};
%get sf and ori curves
%%
figure
for c = 1:length(colordomdum)
    idex = 1;
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
        kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        if c == 4
            kernplot(:,:,:,:) = kernC{1}{p} + kernC{2}{p};
        else
            kernplot(:,:,:,:) = kernC{c}{p};
        end
        
        if f1flag
            %kernplot = 2*mean(phi(TC.tcphase{c}(p,:)).*exp(1i*DM.phasedom*pi/180)); %peak amplitude
            kernplot = 2*abs(mean(kernplot.*phasekern,3)); %peak amplitude
        else
            kernplot = mean(kernplot,3); %average over phase
        end
        
        kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        kernsigplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        if c == 4
            kernsigplot(:,:,:,:) = (kernSigC{1}{p} + kernSigC{2}{p})/2;
        else
            kernsigplot(:,:,:,:) = kernSigC{c}{p};
        end
        kernsigplot = mean(kernsigplot,3); %average over phase
        kernsigplot = reshape(kernsigplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        %kernplot = diff(kernplot,[],3);
        %kernsigplot = diff(kernsigplot,[],3);
 
        %compute optimal time delay
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window   
        maxprof = max(reshape(kerndum(:),length(DM.oridom)*length(DM.sfdom),length(delayWinID(1):delayWinID(2))));
        [dum idma] = max(maxprof); %time delay with max variance
        TC.tauID{c}(p) = idma + delayWinID(1) - 1;
        
        %new sigma after smoothing requires the following
        %Etc2 = kernsigplot.^2 + kernplot.^2;       %get E(x^2) = var(x) + (E(x))^2 
        %kernsigplot = ifftn(fftn(Etc2.^2).*abs(fftn(sfsmoother.^2)));
        %kernsigplot = kernsigplot - (ifftn(fftn(kernplot).*abs(fftn(sfsmoother)))).^2; %This is what my math said
        %kernsigplot = sqrt(kernsigplot);
        %kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(sfsmoother)));
        
        %These are for plotting the time courses in the examples and for
        %establishing significance
        kernTplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));  
        kernsigTplot = ifftn(fftn(kernsigplot).*abs(fftn(tsmoother)));
        
        %This is for making the tuning curves
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother))); %make sure this is done after computing kernsigplot
             
        tcorisfraw = kernplot(:,:,TC.tauID{c}(p));        
        tcorisfrawsig = kernsigplot(:,:,TC.tauID{c}(p)); 
        
        w = phi(mean(tcorisfraw,2))+eps;  %estimate of ori curve 
        %[dum w] = Gaussfit(DM.oridom,w',1); w = w';
        %w = phi(w-max(w(:))/2);
        w = w/sum(w(:));
        tcsf = w'*tcorisfraw; %weighted average over ori
        tcsfsig = w'*tcorisfrawsig; %weighted average over ori
        %tcsf = phi(tcsf);
        
        w = phi(mean(tcorisfraw,1))+eps;  %estimate of spat freq
        %[dum w] = DoGfit(w,DM.sfdom);
        %w = phi(w-max(w(:))/2);
        w = w/sum(w(:));
        tcori = tcorisfraw*w';
        tcorisig = tcorisfrawsig*w';
        %tcori = phi(tcori);

        subdum = DM.sfdom(end)*.01;
        DM.sfdomI = logspace(log10(DM.sfdom(1)),log10(DM.sfdom(end)-subdum),11);     %first interpolate in sp freq   
        if ~isempty(find(isnan(DM.sfdomI)));
            DM.sfdomI = DM.sfdom;
        end
        
        %DM.sfdomI = linspace(DM.sfdom(1),DM.sfdom(end)-.1,14); 
        try
        tcsfI = interp1(DM.sfdom',tcsf',DM.sfdomI','spline')';    
        catch
            'hi'
        end
        
        %[RFenv RFsize] = getRFenvelope(tcsf,DM.sfdom); 
        RFenv = NaN; RFsize = NaN;
        
        if length(DM.oridom)>1
            dI = 5; DM.oridomI = 0:dI:180-dI;
            tcoriI = interp1([DM.oridom 180], [tcori; tcori(1)], [DM.oridomI 180]);  %interpolate in ori
            tcoriI = tcoriI(1:end-1);
            
            tcorisigI = interp1([DM.oridom 180], [tcorisig; tcorisig(1)], [DM.oridomI 180]);  %interpolate in ori
            tcorisigI = tcorisigI(1:end-1);
        else
            DM.oridomI = DM.oridom;
        end        
        
%         if isnan(norm(tcoriI)) || isnan(norm(tcsfI))
%             'hello'
%         end
        
        cothresh = 0.61;

      
 %       [SFparams ffitIsf domI domII pfit Gvaracc repFitQuality] = analyzeLOGtuning(tcsf(:)',DM.sfdom,cothresh);
        %[SFparams ffitIsf domI domII pfit Gvaracc repFitQuality] = analyzeLINtuning(tcsf(:)',DM.sfdom,cothresh);
        
        
         [SFparams ffitIsf domI domII pfit Gvaracc repFitQuality] = analyzeLINtuning_DoG2(tcsf(:)',DM.sfdom,cothresh);
%               
%         [SFparams_L ffitIsf_L domI_L domII_L pfit_L Gvaracc_L repFitQuality_L] = analyzeLOGtuning(tcsf(:)',DM.sfdom,cothresh);
%         if Gvaracc_L > Gvaracc
%             SFparams = SFparams_L ;
%             ffitIsf = ffitIsf_L;
%             domI = domI_L ;
%             domII = domII_L;
%             pfit = pfit_L ;
%             Gvaracc = Gvaracc_L ;
%             repFitQuality = repFitQuality_L;
%         end
        
        
        ffit_sf{p} = ffitIsf;
        DM.sfdomI = domII;
        varacc_sf = Gvaracc;

        TC.sfBWLin{c}(p) = NaN;
        
        %[RFenv RFsize] = getRFenvelope(tcsf,DM.sfdom); 
        
        
        if length(DM.oridomI)>1
            [param ffit_ori{p} varacc_ori sigma] = Gaussfit(DM.oridomI,tcoriI,1);  param(2) = sigma;
            %[param ffit_ori{p} varacc_ori sigma] = Gaussfit_sig(DM.oridomI,tcoriI,tcorisigI,1);  param(2) = sigma;
            %[mu sigma ffit_ori{p} varacc_ori] = CircGaussFit(tcoriI);    param(1) = mu; param(2) = sigma;
        else
            param = NaN*ones(1,4);
            ffit_ori{p} = NaN;
            varacc_ori = 1;
        end

        [idori idsf] = find(tcorisfraw == max(tcorisfraw(:)));
        tcoursema = squeeze(kernTplot(idori,idsf,:));
        tsigcoursema = squeeze(kernsigTplot(idori,idsf,:));
        
        %TC.tcourse{c}(p,:) = tcoursema;

        [idori idsf] = find(tcorisfraw == min(tcorisfraw(:)));
        tcoursemi = squeeze(kernplot(idori,idsf,:));
        tsigcoursemi = squeeze(kernsigplot(idori,idsf,:));
        
        [dum idma] = max(tcoursema);
        %TC.SNR{c}(p) = (tcoursema(idma) - tcoursema(idz))/(tsigcoursema(idma) + tsigcoursema(idz));
        %TC.SNR{c}(p) = (tcoursema(idma))/(tsigcoursema(idma) + sigBlank{p}(idma));

        TC.SNR{c}(p) = (tcoursema(idma) - tcoursemi(idma))/(tsigcoursema(idma) + tsigcoursemi(idma));

        
        TC.tcsfall{c}(p,:) = tcsfI;
        TC.tcoriall{c}(p,:) = tcoriI;
        TC.tcsfall_fit{c}(p,:) = ffit_sf{p};
        TC.tcoriall_fit{c}(p,:) = ffit_ori{p};               
        TC.tcorisfraw{c}(p,:,:) = tcorisfraw;
        
        %I want to keep some values regardless of whether they are lowpass.
        % The others get propagated into map generation and pairwise
        % analysis.
        oriT = 0.0; sfT = 0.0; SNRT = 0.0;
        
        if varacc_ori > oriT & varacc_sf > sfT & TC.SNR{c}(p) > SNRT 
            [TC.flo{c}(p) TC.fhi{c}(p)] = gethcolco(domII,ffitIsf',cothresh);
            TC.sfprefAllpass{c}(p) = SFparams.pref;        
            TC.orisigAllpass{c}(p) = param(2);
            if isnan(TC.flo{c}(p))
                TC.flo{c}(p) = DM.sfdom(1);
            end
        else
            TC.flo{c}(p) = NaN;
            TC.fhi{c}(p) = NaN;
            TC.sfprefAllpass{c}(p) = NaN;
            TC.orisigAllpass{c}(p) = NaN;
        end
            
        
        if varacc_ori > oriT  & varacc_sf > sfT & TC.SNR{c}(p) > SNRT & SFparams.pref>=0.5
            
            %First get all the spatial frequency metrics
            [ma id] = max(ffitIsf);
            
            TC.sfmag{c}(p) = (max(ffitIsf)-min(ffitIsf))/(max(ffitIsf)+min(ffitIsf));
            TC.sfmag{c}(p) = phi(TC.sfmag{c}(p));  %very rarely is it <0
           
            TC.sfpref{c}(p) = SFparams.pref;             
   
            TC.Qfac{c}(p) = TC.sfpref{c}(p)/(TC.fhi{c}(p)-TC.flo{c}(p));
            
            %TC.sfBW{c}(p) = log2(TC.fhi{c}(p)/TC.flo{c}(p)); %bandwidth in octaves               

            TC.sfBWLin{c}(p) = TC.fhi{c}(p)- TC.flo{c}(p); %2sigma bandwidth in cyc/deg from Gaussian fit  
            
            TC.sfBW{c}(p) = log2((TC.sfBWLin{c}(p)/2 + TC.sfpref{c}(p))./TC.sfpref{c}(p))*2;
            
            if TC.sfBW{c}(p) < .2  %high pass cells don't give values at zero, this is a reasonable cutoff                
                TC.sfBW{c}(p) = NaN;                
            end
            
            TC.LPness{c}(p) = (ffitIsf(1))/(max(ffitIsf));           
            
            %TC.LPness{c}(p) = tcsf(1)/max(tcsf);   
            
            %Now get all the orientation metrics
            [TC.OMag{c}(p) TC.OAng{c}(p)] = orifind(tcori,DM.oridom); 
            TC.opref{c}(p) = param(1);
            TC.orisig{c}(p) = param(2);
            %TC.OAng{c}(p) = param(1);
            
            TC.RFsize{c}(p) = RFsize;           
            
            
        else
            %spatial frequency:
            TC.sfpref{c}(p) = NaN;        
            
            TC.sfmag{c}(p) = 0;
            TC.Qfac{c}(p) = NaN;
            TC.sfBW{c}(p) = NaN;
            TC.sfBWsig{c}(p) = NaN;
            TC.sfBWLin{c}(p) = NaN;
            TC.LPness{c}(p) = NaN;
            
            
            %orientation
            TC.OMag{c}(p) = NaN;     
            TC.OAng{c}(p) = NaN; 
            TC.opref{c}(p) = NaN;
            TC.orisig{c}(p) = NaN;
            
            TC.RFsize{c}(p) = NaN;
            
            TC.tcsfall_fit{c}(p,:) = NaN;
            TC.tcoriall_fit{c}(p,:) = NaN;
            
            TC.tcorisfraw{c}(p,:,:) = NaN;

        end                
        
        %This is for the color tuning analyses;  It is redundant for each
        %cell 'c', but it is convenient later
        
        %TC.tccolorall{1}(p,c) = sqrt(max(TC.tcoriall_fit{c}(p,:)) * max(TC.tcsfall{c}(p,:)));      
        TC.tccolorall{1}(p,c) = sqrt(range(TC.tcoriall_fit{c}(p,:)) * range(TC.tcsfall{c}(p,:)));  
        %TC.tccolorall{1}(p,c) = var(kernplot(:))*sign(TC.opref{c}(p));
        
        
%         semilogx(DM.sfdom,tcsf,['.' colorid{c} '-']), hold on
%         plot(DM.sfdomI,ffitIsf,'k')
%         title(num2str(TC.sfBWLin{c}(p)))
% 
%         axis tight
        
        imagesc(tcorisfraw)
        axis off
        drawnow

                %errorbar(DM.oridom,tcori,tcorisig,['.' colorid{c} '-']), hold on
%                 plot(DM.oridom,tcori,['.' colorid{c} '-']), hold on
%                 plot(DM.oridomI,ffit_ori{p},'k')
%                 axis tight

        %
                 title([num2str(round(TC.orisig{c}(p))) ' ' num2str(TC.sfpref{c}(p))])

%        imagesc(tcorisfraw)
%         hold on
%         contour(ffitorisf,(max(ffitorisf(:))+min(ffitorisf(:)))/2,'k')
         %title(num2str(TC.sfpref{c}(p)))
         
         %Get stuff to plot the examples
         if ~isempty(find(p == idExamp))
            kern_examp{c}{idex} = tcorisfraw;   
            
            [idori idsf] = find(tcorisfraw == max(tcorisfraw(:)));
            tcoursema_examp{c}{idex} = squeeze(kernTplot(idori,idsf,:));            
            tsigcoursema_examp{c}{idex} = squeeze(kernsigTplot(idori,idsf,:));
            
            [idori idsf] = find(tcorisfraw == min(tcorisfraw(:)));
            tcoursemi_examp{c}{idex} = squeeze(kernTplot(idori,idsf,:));            
            tsigcoursemi_examp{c}{idex} = squeeze(kernsigTplot(idori,idsf,:));
            
            tcoriEx{c}{idex} = tcori;
            tcsfEx{c}{idex} = tcsf;
            
            ffitoriEx{c}{idex} = ffit_ori{p};
            ffitsfEx{c}{idex} = ffitIsf;
            
            hh = makeTemporalfilter;
            acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;
            
            mu = mean(cellS.cellMat{1}(idex+1,:));
            y = processTcourse(cellS.cellMat{1}(idex+1,:),hh,1,acqPeriod)+mu; 
            id = find(y<prctile(y,50) & y>prctile(y,0));
            Fo = median(y(id));
            tcourse{c}{idex} = (y-Fo)/Fo; 
            
            %This is for Eyal's paper
            oridynamicsEx{c}{idex} = squeeze(mean(kernTplot(:,3:end-1,:),2));
            
            idex = idex+1;
         end
        
    end

end



%% Compute F1/F0

sig = 50;
dom = DM.taudom-mean(DM.taudom);
psmooth = exp(-dom.^2/(2*sig^2));
psmooth = psmooth/sum(psmooth);
psmooth = ones(length(DM.phasedom),1)*psmooth;
psmooth = abs(fft(psmooth,[],2));

phasedomI = 0:359;
%kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));
figure
for c = 1:length(colordomdum)
    idex = 1;
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
        
        %if its a color stimulus, use the first and second color (L&M) to
        %get preferred ori
        
        if length(colordomdum) > 1
            muori = angle(exp(1i*TC.OAng{1}(p)*pi/180) + exp(1i*TC.OAng{2}(p)*pi/180))*180/pi;
            if muori<0
                muori = muori+180;
            end
        else
            %muori = nanmean(exp(1i*2*TC.OAng{c}(:)*pi/180));
            %muori = angle(muori)*180/pi/2 ; %Use one ori for all cells
            
            
            muori = TC.OAng{c}(p);  %Best ori for this cell

        end     
        
        
        base = prctile(kernC{c}{p}(:),30);
        if c == 4
            kdum = (kernC{1}{p} + kernC{2}{p})/2; 
        else
            kdum = kernC{c}{p}; 
        end
        
        domdum = DM.oridom;
        if length(DM.oridom) >= 12
            kdum = (kdum(1:2:end,:,:,:) + kdum(2:2:end,:,:,:))/2;
            domdum = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
        end        
        
        [dum oriid] = min(abs(muori-domdum));
        if length(colordomdum) > 1
            musf = sqrt(TC.sfpref{1}(p).*TC.sfpref{2}(p));
        else 
            musf = TC.sfpref{1}(p); %Use best sf for this cell
            
            %musf = 10^nanmean(log10(TC.sfpref{1}(:))); %Use one sf for all cells
            
        end
        [dum sfid] = min(abs(musf-DM.sfdom));

        kernplot = squeeze(kdum(oriid,sfid,:,:));  %phase and time
        if size(kernplot,2) == 1
            kernplot = kernplot';
        end
        kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2); %smooth over time
        
        TC.tcphase{c}(p,:) = squeeze(kernplot(:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves

        f1 = 2*mean(phi(TC.tcphase{c}(p,:)).*exp(1i*DM.phasedom*pi/180)); %peak amplitude
        f0 = mean(phi(TC.tcphase{c}(p,:)));

        TC.phase{c}(p) = angle(f1)*180/pi;       

        TC.F1F0{c}(p) = 2*abs(f1)/f0;  %peak-to-peak amplitude, over mean
        
    
        phasefitI = abs(f1)*cos(phasedomI*pi/180 - TC.phase{c}(p)*pi/180) + f0;
        phasefit = abs(f1)*cos(DM.phasedom*pi/180 - TC.phase{c}(p)*pi/180) + f0;
        
        %F1 as measured by Nishimoto
         %dum = phi(TC.tcphase{c}(p,:));
         %f1dum = abs(dum(1)-dum(3)) + abs(dum(2)-dum(4));
         %TC.F1F0{c}(p) = f1dum/mean(dum);
        
        
        varacc = (var(TC.tcphase{c}(p,:))-var(TC.tcphase{c}(p,:)-phasefit))/var(TC.tcphase{c}(p,:));
        
        if varacc < .7
        %if TC.F1F0{c}(p) < .8
            TC.phase{c}(p) = NaN;
            %TC.F1F0{c}(p) = NaN;
        end
        if isnan(TC.opref{c})
        %if TC.F1F0{c}(p) < .8
            TC.phase{c}(p) = NaN;
            TC.F1F0{c}(p) = NaN;
        end
        

        %The Nishimoto et al version 
%         f1 = abs(TC.tcphase{c}(p,1)-TC.tcphase{c}(p,3)) + abs(TC.tcphase{c}(p,2)-TC.tcphase{c}(p,4));
%         f0 = mean(TC.tcphase{c}(p,:));
%         TC.F1F0{c}(p) = abs(f1)/f0;

        plot([DM.phasedom 360],[TC.tcphase{c}(p,:) TC.tcphase{c}(p,1)] ,['.-' colorid{c}]), hold on
        plot(phasedomI,phasefitI,'k'), 
        hold on, plot([0 360],[0 0],'--')
        axis tight

        axis off
        %title(num2str(TC.phase{c}(p)))
        title(num2str(round(TC.F1F0{c}(p)*10)/10))
        if isnan(TC.F1F0{c}(p))
           title('') 
        end

        %Get stuff to plot the examples
        if ~isempty(find(p == idExamp))

            tcphaseEx{c}{idex} = TC.tcphase{c}(p,:);
           
            ffitphaseEx{c}{idex} = phasefit;
            
            ffitphaseExI{c}{idex} = phasefitI;

            idex = idex+1;
        end

    end

end


for c = 1:length(colordomdum)
    id = find(isnan(TC.OMag{c}));
    TC.F1F0{c}(id) = NaN;
end


%id = find(TC.F1F0{c}<0);
%TC.F1F0{c}(id) = 0;
%id = find(TC.F1F0{c}>2);
%TC.F1F0{c}(id) = 2;
figure,hist(TC.F1F0{1},linspace(0,3,20))

%% Get phase bias for each ori sf combination

% sig = 50;
% dom = DM.taudom-mean(DM.taudom);
% psmooth = exp(-dom.^2/(2*sig^2));
% psmooth = psmooth/sum(psmooth);
% psmooth = ones(length(DM.phasedom),1)*psmooth;
% psmooth = abs(fft(psmooth,[],2));
% 
% phasedomI = 0:359;
% phistdom = [-180:45:180];
% 
% sigmat = zeros(length(DM.sfdom),length(DM.oridom));
% clear countP phaseSelectivityBoot
% c = 1;
% k = 1;
% figure
% for sfid = 1:length(DM.sfdom)
%     for oriid = 1:length(DM.oridom)        
%         
%         for p = 1:MK.Ncell
%             
%             %if its a color stimulus, use the first and second color (L&M) to
%             %get preferred ori
%             
%             muori = DM.oridom(oriid);
%             
%             base = prctile(kernC{c}{p}(:),30);
%             
%             kdum = kernC{c}{p};
%             
% %             domdum = DM.oridom;
% %             if length(DM.oridom) >= 6
% %                 kdum = (kdum(1:2:end,:,:,:) + kdum(2:2:end,:,:,:))/2;
% %                 domdum = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
% %             end
% %             
% %             [dum oriid] = min(abs(muori-domdum));
%             
%             musf = DM.sfdom(sfid);
%             
%             %[dum sfid] = min(abs(musf-DM.sfdom));
%             
%             kernplot = squeeze(kdum(oriid,sfid,:,:));  %phase and time
%             if size(kernplot,2) == 1
%                 kernplot = kernplot';
%             end
%             kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2); %smooth over time
%             
%             
%             tcdum = squeeze(kernplot(:,TC.tauID{c}(p)))';  %use same time delay as ori/sf curves
%             
%             f1 = 2*mean(phi(tcdum).*exp(1i*DM.phasedom*pi/180)); %peak amplitude
%             f0 = mean(phi(tcdum));
%             
%             phasepref{oriid,sfid}(p) = angle(f1)*180/pi;
%             phaseamp{oriid,sfid}(p) = abs(f1);
%             
%             phasefitI = abs(f1)*cos(phasedomI*pi/180 - phasepref{oriid,sfid}(p)*pi/180) + f0;
%             phasefit = abs(f1)*cos(DM.phasedom*pi/180 - phasepref{oriid,sfid}(p)*pi/180) + f0;
%             
%             phasefitAll{oriid,sfid}(p,:) = phasefitI;
%             tcphaseAll{oriid,sfid}(p,:) = tcdum;
%             
%             odiff = abs(oridiff(TC.opref{1}(p)*pi/180,DM.oridom(oriid)*pi/180))*180/pi;
%             odiffnorm = odiff/TC.orisig{1}(p)
%             
%             %F1 as measured by Nishimoto
%             %         dum = phi(TC.tcphase{c}(p,:));
%             %         f1dum = abs(dum(1)-dum(3)) + abs(dum(2)-dum(4));
%             %         f1fodum = f1dum/mean(dum);
%             
%             varacc = (var(tcdum)-var(tcdum-phasefit))/var(tcdum);
%             
%             f1f0 = 2*phaseamp{oriid,sfid}(p)/f0;
%             
%             if varacc > .7 & f1f0 > .5 % & odiffnorm<3
%                 
%                 yieldMat{oriid,sfid}(p) = 1;
%             else
%                 
%                 %if TC.F1F0{c}(p) < .8
%                 phasepref{oriid,sfid}(p) = NaN;
%                 phaseamp{oriid,sfid}(p) = NaN;
%                 
%                 yieldMat{oriid,sfid}(p) = 0;
% 
%             end
%             
%             %The Nishimoto et al version
%             %         f1 = abs(TC.tcphase{c}(p,1)-TC.tcphase{c}(p,3)) + abs(TC.tcphase{c}(p,2)-TC.tcphase{c}(p,4));
%             %         f0 = mean(TC.tcphase{c}(p,:));
%             %         TC.F1F0{c}(p) = abs(f1)/f0;
%        
%             
%         end
%         
%         subplot(length(DM.sfdom),length(DM.oridom),k)
%         k = k+1;
%         
%         phaseprefdum = phasepref{oriid,sfid};
%         %muang = angle(nanmean(exp(1i*phaseprefdum*pi/180)))*180/pi;
%         %phaseprefdum = angle(exp(1i*(phaseprefdum-muang)*pi/180))*180/pi;
%         
%         
%         
%          h = hist(phaseprefdum,phistdom);
%         %polar(phistdom*pi/180,h,'o-k')
%         bar(phistdom,h,'FaceColor',[1 1 1])
%         
%         %polar(phistdom,h)
%         %axis off
%         ylimheight = 25;
%         ylim([0 ylimheight])
%         axis off
%         
%          Res = sum(h.*exp(1i*phistdom*pi/180));
%          Resang = angle(Res)*180/pi;
%          phaseSelectivity(sfid,oriid) = abs(Res)/sum(h);
%          hold on
%          plot([Resang Resang],[0 phaseSelectivity(sfid,oriid)*ylimheight],'r','LineWidth',2)
%          
%          %Does each histogram yield a phase selectivity value that is
%          %statististically unique from a random sampling from a uniform
%          %distribution with the same number of data points?
%          Niter = 1000;
%          for iter = 1:Niter  
%              
%             phaseprefshuff = 360*rand(1,length(phasepref{oriid,sfid}))-180;
%             phaseprefshuff(find(isnan(phasepref{oriid,sfid}))) = NaN;
% 
%            
%             h = hist(phaseprefshuff,phistdom);
%             
%             Res = sum(h.*exp(1i*phistdom*pi/180));           
% 
%             phaseSelectivityBoot(sfid,oriid,iter) = abs(Res)/sum(h);
%             
%             countP(iter) = (sign(phaseSelectivity(sfid,oriid) - phaseSelectivityBoot(sfid,oriid,iter))+1)/2;
%          end
%          if sum(countP)/Niter>.99
%              title('*')
%              sigmat(sfid,oriid) = 1;
%          end
%          
%         %plot(phasedomI,phasefitAll{oriid,sfid}','k'), axis tight
%         
%         
%     end
% end
% 
% figure,imagesc(DM.oridom,DM.sfdom,phaseSelectivity,[0 .7])
% xlabel('ori'), ylabel('sf')
% colorbar
% 
% sfrange = 4:6;  %Should not include low sfs
% phaseSelBootDum = phaseSelectivityBoot(sfrange,:,:);
% phaseSelDum = phaseSelectivity(sfrange,:);
% 
% muS = mean(phaseSelBootDum,3);
% sS = std(phaseSelBootDum,[],3);
% 
% figure,
% subplot(2,1,1),hist((phaseSelDum(:)-muS(:))./sS(:),8)
% xlabel('(Resultant - E(Resultant))/std(Resultant)')
% subplot(2,1,2),imagesc(DM.oridom(:),DM.sfdom(sfrange),(phaseSelDum-muS)./sS)
% xlabel('ori'), ylabel('sf')
% colorbar
% 
% [h p] = ttest((phaseSelDum(:)-muS(:))./sS(:))
% 
% figure,scatter(phasepref{2,1},phasepref{2,2})


%%
% figure
% k = 1;
% 
% 
% 
% for sfid = 1:length(DM.sfdom)
%     
%     RetEnvAll{sfid} = [];
%     phidumAll{sfid} = [];
%     for oriid = 1:length(DM.oridom)               
%         
%         %%%%Get envelope location from low spatial frequencies at this
%         %%%%orientation
%         sfrefA = 1;
%         Res = (nanmean(exp(1i*phasepref{oriid,sfrefA}*pi/180)));
%         muphi = angle(Res);
%         phidumA = angle(exp(1i*(phasepref{oriid,sfrefA}*pi/180-muphi)))*180/pi;
%         
%         sfrefB = 2;
%         Res = (nanmean(exp(1i*phasepref{oriid,sfrefB}*pi/180)));
%         muphi = angle(Res);
%         phidumB = angle(exp(1i*(phasepref{oriid,sfrefB}*pi/180-muphi)))*180/pi;
%         
%         RetEnv = (phidumA/DM.sfdom(sfrefA)/360 + phidumB/DM.sfdom(sfrefB)/360)/2;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         Res = (nanmean(exp(1i*phasepref{oriid,sfid}*pi/180)));
%         muphi = angle(Res);
%         phidum = angle(exp(1i*(phasepref{oriid,sfid}*pi/180-muphi)))*180/pi;
%         
%         subplot(length(DM.sfdom),length(DM.oridom),k)
%         scatter(RetEnv,phidum/DM.sfdom(sfid)/360,'.k')
%         %scatter(phidumLow,phidum), xlim([-180 180]), ylim([-180 180])
%         xlim([-1 1]), 
%         
%         speriod = round(1/DM.sfdom(sfid)*100)/100;
%         
%         ylim([-speriod/2 speriod/2])
%         ylim([-1 1]), 
%         
%         hold on,
%         plot([-2 2],[speriod/2 speriod/2])
%         hold on,
%         plot([-2 2],[-speriod/2 -speriod/2])
%         
%         hold on
%         plot([-2 2],[-2 2])
%         
%         set(gca,'Ytick',[-speriod/2 speriod/2])
%         
%         axis square
%         
%          
%         
%         k = k+1;        
%         
%         RetEnvAll{sfid} = [RetEnvAll{sfid}; RetEnv(:)];
%         phidumAll{sfid} = [phidumAll{sfid}; phidum(:)/DM.sfdom(sfid)/360];
%         
%     end
% end

%%

% figure
% for sfid = 1:length(DM.sfdom)
%     subplot(length(DM.sfdom),1,sfid)
%     scatter(RetEnvAll{sfid},phidumAll{sfid},'.k')
%     
%     hold on
%     plot([-2 2],[-2 2])
%     ylim([-1 1]),
%     xlim([-1 1]),
%     
%             speriod = round(1/DM.sfdom(sfid)*100)/100;
%         
%         ylim([-speriod/2 speriod/2])
% 
%         
%         hold on,
%         plot([-2 2],[speriod/2 speriod/2])
%         hold on,
%         plot([-2 2],[-speriod/2 -speriod/2])
%     
% end

%% Plot phase map for each ori/sf combo
% [xmicperpix ymicperpix] = getImResolution;
% xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
% ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
% 
% 
% 
% k = 1;
% figure
% for sfid = 1:length(DM.sfdom)
%     for oriid = 1:length(DM.oridom)
%         
%         phaseprefdum = phasepref{oriid,sfid};
%         muang = angle(nanmean(exp(1i*phaseprefdum*pi/180)))*180/pi;
%         phaseprefdum = angle(exp(1i*(phaseprefdum-muang)*pi/180))*180/pi;
%         
%         phaseprefdum = phasepref{oriid,sfid};
%         
%         s = subplot(length(DM.sfdom),length(DM.oridom),k);
%         plotDotMap(phaseprefdum,'hsv',MK.masklabel,MK.celldom(MK.nID),[-180 180]);
%         phasedom = {-180 ,-90 ,0 ,90, 180};
%         
%         
%         %         plotDotMap(phaseprefdum/DM.sfdom(sfid)/360,'hsv',MK.masklabel,MK.celldom(MK.nID),[-180 180]/DM.sfdom(sfid)/360);
%         %         degTicks = [-180 -90 0 90 180]/360/DM.sfdom(sfid);
%         %         degTicks = round(degTicks*100)/100;
%         %         phasedom = {degTicks(1) ,degTicks(2) ,degTicks(3) ,degTicks(4), degTicks(5)};  %units of degrees
%         
%         if sfid == 1 & oriid == 1
%             iddom = linspace(0,1,length(phasedom));
%             colorbar('YTick',iddom,'YTickLabel',phasedom)
%         end
%         %axis off
%         set(gca,'XTick', [], 'YTick',[])
%         colormap(s,hsv)
%         %title('phase preference')
%         if sigmat(sfid,oriid)
%             title('*')
%         end
%         
%         k = k+1;
%     end
% end
% 
% 
% 

%% Bootstrap; This shuffles each tuning curve
% c = 1;
% k = 1;
% phaseSelectivityBoot = 0;
% figure
% for iter = 1:100
%     iter
%     for sfid = 1:length(DM.sfdom)
%         for oriid = 1:length(DM.oridom)
%             
%             for p = 1:MK.Ncell
%                 
%                 %if its a color stimulus, use the first and second color (L&M) to
%                 %get preferred ori
%                 
%                 muori = DM.oridom(oriid);
%                 
%                 base = prctile(kernC{c}{p}(:),30);
%                 
%                 kdum = kernC{c}{p};
%                 
%                 %             domdum = DM.oridom;
%                 %             if length(DM.oridom) >= 6
%                 %                 kdum = (kdum(1:2:end,:,:,:) + kdum(2:2:end,:,:,:))/2;
%                 %                 domdum = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
%                 %             end
%                 %
%                 %             [dum oriid] = min(abs(muori-domdum));
%                 
%                 musf = DM.sfdom(sfid);
%                 
%                 %[dum sfid] = min(abs(musf-DM.sfdom));
%                 
%                 kernplot = squeeze(kdum(oriid,sfid,:,:));  %phase and time
%                 if size(kernplot,2) == 1
%                     kernplot = kernplot';
%                 end
%                 kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2); %smooth over time
%                 
%                 tcdum = squeeze(kernplot(:,TC.tauID{c}(p)))';  %use same time delay as ori/sf curves
%                 
%                 %%%%%Shuffle%%%%%%
%                 [dum id] = sort(rand(1,length(tcdum)));
%                 tcdum = tcdum(id);
%                 %%%%%%%%%%%%%%%%%
%                 
%                 f1 = 2*mean(phi(tcdum).*exp(1i*DM.phasedom*pi/180)); %peak amplitude
%                 f0 = mean(phi(tcdum));
%                 
%                 phasepref{oriid,sfid}(p) = angle(f1)*180/pi;
%                 phaseamp{oriid,sfid}(p) = abs(f1);
%                 
%                 phasefitI = abs(f1)*cos(phasedomI*pi/180 - phasepref{oriid,sfid}(p)*pi/180) + f0;
%                 phasefit = abs(f1)*cos(DM.phasedom*pi/180 - phasepref{oriid,sfid}(p)*pi/180) + f0;
%                 
%                 phasefitAll{oriid,sfid}(p,:) = phasefitI;
%                 tcphaseAll{oriid,sfid}(p,:) = tcdum;
%                 
%                 %F1 as measured by Nishimoto
%                 %         dum = phi(TC.tcphase{c}(p,:));
%                 %         f1dum = abs(dum(1)-dum(3)) + abs(dum(2)-dum(4));
%                 %         f1fodum = f1dum/mean(dum);
%                 
%                 varacc = (var(tcdum)-var(tcdum-phasefit))/var(tcdum);
%                 
%                 if ~yieldMat{oriid,sfid}(p);
%                     %if TC.F1F0{c}(p) < .8
%                     phasepref{oriid,sfid}(p) = NaN;
%                     phaseamp{oriid,sfid}(p) = NaN;
%                 end
%                 
%                 %The Nishimoto et al version
%                 %         f1 = abs(TC.tcphase{c}(p,1)-TC.tcphase{c}(p,3)) + abs(TC.tcphase{c}(p,2)-TC.tcphase{c}(p,4));
%                 %         f0 = mean(TC.tcphase{c}(p,:));
%                 %         TC.F1F0{c}(p) = abs(f1)/f0;
%                 
%                 
%             end
%             
%             %subplot(length(DM.sfdom),length(DM.oridom),k)
%             k = k+1;
%             
%             h = hist(phasepref{oriid,sfid},phistdom);
%             %polar(phistdom*pi/180,h,'o-k')
%             %bar(phistdom,h)
%             %axis off
%             %ylim([0 15])
%             
%             Res = sum(h.*exp(1i*phistdom*pi/180));
%             
%             
%             %plot(phasedomI,phasefitAll{oriid,sfid}','k'), axis tight
%             phaseSelectivityBoot(sfid,oriid,iter) = abs(Res)/sum(h);
%             
%         end
%     end
% end
% 
% clear countP
% for iter = 1:100
%     for sfid = 1:length(DM.sfdom)
%         for oriid = 1:length(DM.oridom)
%             countP(sfid,oriid,iter) = (sign(phaseSelectivity(sfid,oriid) - phaseSelectivityBoot(sfid,oriid,iter))+1)/2;
%         end
%     end
% end
% 
% pMat = sum(countP,3)/size(countP,3);
% %figure,imagesc(pMat)
% 
% sigMat = zeros(size(pMat));
% sigMat(find(pMat>.99)) = 1;
% 
% %%
% figure,imagesc(sigMat), colormap gray, title('significant clustering')
% 
% figure
% subplot(3,1,1),imagesc(DM.oridom,DM.sfdom,phaseSelectivity,[0 1])
% xlabel('ori'), ylabel('sf')
% subplot(3,1,2),
% imagesc(DM.oridom,DM.sfdom,mean(phaseSelectivityBoot,3),[0 1])
% title('Bootstrap phase selecitivity (mean)')
% subplot(3,1,3),
% imagesc(DM.oridom,DM.sfdom,std(phaseSelectivityBoot,[],3),[0 1])
% title('Bootstrap phase selecitivity (std)')
% 
% xlabel('ori'), ylabel('sf')
% colorbar
% %% Get a single stat across all ori and sfs
% sfrange = 4:5;  %Should not 
% phaseSelBootDum = phaseSelectivityBoot(sfrange,:,:);
% phaseSelDum = phaseSelectivity(sfrange,:);
% 
% muS = mean(phaseSelBootDum,3);
% sS = std(phaseSelBootDum,[],3);
% 
% figure,
% subplot(2,1,1),hist((phaseSelDum(:)-muS(:))./sS(:),8)
% xlabel('(Resultant - E(Resultant))/std(Resultant)')
% subplot(2,1,2),imagesc(DM.oridom(:),DM.sfdom(sfrange),(phaseSelDum-muS)./sS)
% xlabel('ori'), ylabel('sf')
% colorbar
% 
% [h p] = ttest((phaseSelDum(:)-muS(:))./sS(:))

%%
% 
% figure,
% subplot(1,2,1), 
% phistdom = [-180:45:180];
% h = hist(TC.phase{1},phistdom);
% polar(phistdom*pi/180,h,'o-k')
% 
% subplot(1,2,2), 
% 
% tcphasenorm = TC.tcphase{c}';
% ma = ones(size(tcphasenorm,1),1)*max((tcphasenorm));
% mi = ones(size(tcphasenorm,1),1)*min((tcphasenorm));
% tcphasenorm = (tcphasenorm-mi)./(ma-mi);
% plot(DM.phasedom,tcphasenorm,'.-k'), 
% %hold on
% 
% o = 2
% s = 5;
% tcphasenorm = phasefitAll{o,s}';
% ma = ones(size(tcphasenorm,1),1)*max((tcphasenorm));
% mi = ones(size(tcphasenorm,1),1)*min((tcphasenorm));
% tcphasenorm = (tcphasenorm-mi)./(ma-mi);
% 
% plot(phasedomI,tcphasenorm','r'), axis tight
% 
% for c = 1:length(colordomdum)
%     id = find(isnan(TC.OMag{c}));
%     TC.F1F0{c}(id) = NaN;
% end
% 
% 
% %id = find(TC.F1F0{c}<0);
% %TC.F1F0{c}(id) = 0;
% %id = find(TC.F1F0{c}>2);
% %TC.F1F0{c}(id) = 2;
% figure,hist(TC.F1F0{1},linspace(0,3,20))
% 


%% Reconstruct RF
% 
% sig = 80;
% dom = DM.taudom-mean(DM.taudom);
% psmooth = exp(-dom.^2/(2*sig^2));
% psmooth = psmooth/sum(psmooth);
% psmooth = abs(fft(psmooth,[],2));
% 
% 
% 
% sigmat = zeros(length(DM.sfdom),length(DM.oridom));
% clear countP phaseSelectivityBoot
% c = 1;
% k = 1;
% figure
% 
% 
% for p = 1:MK.Ncell
%     
%     
%     base = prctile(kernC{c}{p}(:),30);
%     
%     kdum = kernC{c}{p};
%    
%    
%     
%     for ori = 1:length(DM.oridom)
%         for sf = 1:length(DM.sfdom)
%             for ph = 1:length(DM.phasedom)
%                
%                 tdum = squeeze(kdum(ori,sf,ph,:))';  %get time trace
%                 tdum = ifft(fft(tdum,[],2).*psmooth,[],2); %Smooth it
%                 kslice{p}(ori,sf,ph)  = tdum(TC.tauID{c}(p));  %take single delay
%                 
%             end
%         end
%     end
%     
% 
%     
% 
%     
%     subplot(10,10,p),
%     imagesc(mean(kslice{p},3))
%     
%     
% end
% 
% RF = GMakeRF2(kslice);
% RFexamp = RF(idExamp);
% 

%%
%Get parameters that combine across colors
if length(DM.colordom)>1
    
    for p = 1:length(TC.tccolorall{1}(:,1))
        TC.tccolorall{1}(p,:) = TC.tccolorall{1}(p,:)/max(TC.tccolorall{1}(p,:));

        %lumpref(p) = log((colorvec(p,2)/colorvec(p,3))); %log((L-M)/(L+M)) OR log(M/S)       
        
        if strcmp(getparam('colorspace'),'LMS')
            id1 = 1; id2 = 2;
        elseif strcmp(getparam('colorspace'),'DKL')
            id1 = 1; id2 = 3;
        end
        
        TC.lumpref{1}(p) = log2((TC.tccolorall{1}(p,id1)/TC.tccolorall{1}(p,id2))); %log(S/(L-M)) OR log(L/M)
        %TC.lumpref{1}(p) = (TC.tccolorall{1}(p,1)-TC.tccolorall{1}(p,2))/(TC.tccolorall{1}(p,1)+TC.tccolorall{1}(p,2)); %L-M/L+M

        TC.LMphaseDiff{1}(p) = angle(exp(1i*pi/180*(TC.phase{id1}(p)-TC.phase{id2}(p))))*180/pi;  %(angle(L) - angle(M))
        
        %LMproj is a "normalized" measure of the opponency.  Its
        %cos(phasedifferenc), then multiplied by the modulation depths:
        %sqrt(F1_L/F0_L * F1_L/F0_L)/2... so it will be -1 for oppenent
        %cells with large modulation depth, +1 for non-opponent cells with
        %large modulation, and ~zero for other cases

        TC.LMproj{1}(p) = sqrt(TC.F1F0{1}(p)*TC.F1F0{2}(p))/2*cos(TC.LMphaseDiff{1}(p)); 
        
        TC.LMphaseDiff{1}(p) = abs(TC.LMphaseDiff{1}(p));

    end
       
    %This is dumb, but convenient later
    TC.tccolorall{1} = TC.tccolorall{1}(:,1:3);
    TC.tccolorall{2} = TC.tccolorall{1};
    TC.tccolorall{3} = TC.tccolorall{1};
    TC.tccolorall{4} = TC.tccolorall{1};
    
    TC.lumpref{1} = TC.lumpref{1};
    TC.lumpref{2} = TC.lumpref{1};
    TC.lumpref{3} = TC.lumpref{1};
    TC.lumpref{4} = TC.lumpref{1};
    
    TC.LMphaseDiff{1} = TC.LMphaseDiff{1};
    TC.LMphaseDiff{2} = TC.LMphaseDiff{1};
    TC.LMphaseDiff{3} = TC.LMphaseDiff{1};
    TC.LMphaseDiff{4} = TC.LMphaseDiff{1};
    
    TC.LMproj{1} = TC.LMproj{1};
    TC.LMproj{2} = TC.LMproj{1};
    TC.LMproj{3} = TC.LMproj{1};
    TC.LMproj{4} = TC.LMproj{1};
    
end


%% Compute the spatial phase for a single ORI and SF
% clear oriid sfid
% for c = 1:length(colordomdum)
%     popOpref = angle(nansum(exp(1i*TC.OAng{1}*2*pi/180)))/2;
%     if popOpref<0;
%         popOpref = popOpref + 180;
%     end
% 
%     id = find(~isnan(TC.sfpref{c}));
%     popSpref = geomean(TC.sfpref{c}(id));
% 
%     [dum oriid{c}] = min(abs(popOpref-DM.oridom));
%     [dum sfid{c}] = min(abs(popSpref-DM.sfdom));
% end
% 
% for p = 1:MK.Ncell
% 
%     for c = 1:length(colordomdum)
%         
%         kernplot = squeeze(kernC{c}{p}(oriid{c},sfid{c},:,:));  %phase and time
%         kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2);
%         
%         tcphase = squeeze(kernplot(:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves
%         
%         tcphase = phi(tcphase);        
%         f1 = sum(tcphase.*exp(1i*DM.phasedom'*pi/180));
%         
%         TC.phase{c}(p) = angle(f1)*180/pi;
%         
%         %The Nishimoto et al version
% %         f1 = abs(tcphase(1)-tcphase(3)) + abs(tcphase(2)-tcphase(4));
% %         f0 = sum(tcphase);
% %         TC.F1F0{c}(p) = 2*f1/f0;
%         
%     end
% 
% end
% 
% for c = 1:length(colordomdum)
%     id = find(isnan(TC.OMag{c}));
%     TC.phase{c}(id) = NaN;
% end

%% Compute 1D spatial profile

% sig = 50;
% dom = DM.taudom-mean(DM.taudom);
% psmooth = exp(-dom.^2/(2*sig^2));
% psmooth = ones(4,1)*psmooth;
% psmooth = abs(fft(psmooth,[],2));
% for i = 1:length(DM.sfdom)
%     smoother(i,:,:) = psmooth;
% end
% 
% figure
% for p = 1:MK.Ncell
%     subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
% 
%     for c = 1:length(colordomdum)
%         
%         [dum oriid] = min(abs(TC.OAng{c}(p)-DM.oridom));
%         
%         kernplot = squeeze(kernC{c}{p}(oriid,:,:,:));  %spatial frequency, phase, and time
%         kernplot = ifft(fft(kernplot,[],3).*smoother,[],3);
%         
%         LineF = squeeze(kernplot(:,:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves
%         
%         
%         
%         plot(DM.phasedom,tcphase,['.' colorid{c} '-']), hold on
%         
%         title(num2str(TC.F1F0{c}(p)))
%         
%     end
% 
% end

%% Plot examples
for i = 1:length(DM.sfdomI)
    sfdomcell{i} = num2str(num2str([round(DM.sfdomI(i)*10)/10]));
end

if ~isempty(idExamp)
    for c = 1:length(DM.colordom)
        figure
        Nex = length(kern_examp{1});
        for i = 1:Nex
            subplot(Nex*2,3,i*6-5)
            imagesc(kern_examp{c}{i}), colormap jet
            %xl = str2num(cell2mat(get(gca,'XtickLabel')));
            %set(gca,'XtickLabel',round(DM.sfdom(xl)*10)/10)
            set(gca,'Xtick',[1 length(DM.sfdom)])
            set(gca,'XtickLabel',round([DM.sfdom(1) DM.sfdom(end)]*10)/10,'TickDir','out')
            set(gca,'Ytick',[1 length(DM.oridom)])
            set(gca,'YtickLabel',round([DM.oridom(1) DM.oridom(end)]*10)/10,'TickDir','out')
            xlabel('sf'), ylabel('ori')
            axis image
            
            subplot(Nex*2,3,i*6-4)

            fill([DM.taudom fliplr(DM.taudom)],[tcoursema_examp{c}{i}-tsigcoursema_examp{c}{i}; flipud(tcoursema_examp{c}{i}+tsigcoursema_examp{c}{i})]',[.0 .0 1])
            hold on
            plot(DM.taudom,tcoursema_examp{c}{i},'k'),xlabel('ms')
            %ylim([min(tcoursema_examp{c}{i}-tsigcoursema_examp{c}{i})-.1 max(tcoursema_examp{c}{i}+tsigcoursema_examp{c}{i})+.1])

            fill([DM.taudom fliplr(DM.taudom)],[tcoursemi_examp{c}{i}-tsigcoursemi_examp{c}{i}; flipud(tcoursemi_examp{c}{i}+tsigcoursemi_examp{c}{i})]',[1 .0 .0])
            hold on
            plot(DM.taudom,tcoursemi_examp{c}{i},'k'),xlabel('ms')
            xlim([DM.taudom(1) DM.taudom(end)])
            axis tight

            subplot(Nex*2,3,i*6-2)
            plot(DM.oridom,tcoriEx{c}{i},'.k'),xlabel('ori')
            hold on, plot(DM.oridomI,ffitoriEx{c}{i},'r')
            set(gca,'Xtick',[0 90 180]), xlim([-5 180])
            %ylim([min(ffitoriEx{c}{i})-.1 max(ffitoriEx{c}{i})+.1])
            axis tight
            yL = get(gca,'ylim');
            if yL(1)>0
                yL(1) = 0;
            end
            yL(2) = yL(2)+.1;
            ylim(yL)

            subplot(Nex*2,3,i*6-1)
            plot(log2(DM.sfdom),tcsfEx{c}{i},'.k'), xlabel('sf')
            hold on, plot(log2(DM.sfdomI),ffitsfEx{c}{i},'r'), %xlim([DM.sfdomI(1) DM.sfdomI(end)])
            %ylim([min(ffitsfEx{c}{i})-.1 max(ffitsfEx{c}{i})+.1])
            set(gca,'XTick',log2(DM.sfdom(1:2:end)))
            set(gca,'XTickLabel',round(DM.sfdom(1:2:end)*10)/10)
            xlim([log2(DM.sfdom(1)-.1) log(DM.sfdom(end))+1.5])
            axis tight
            yL = get(gca,'ylim');
            if yL(1)>0
                yL(1) = 0;
            end
            yL(2) = yL(2)+.1;
            ylim(yL)
            
            subplot(Nex*2,3,i*6)
            plot([DM.phasedom 360],[tcphaseEx{c}{i} tcphaseEx{c}{i}(1)],'.k'),xlabel('phase')
            hold on, plot(phasedomI,ffitphaseExI{c}{i},'r')
            set(gca,'Xtick',[0 180 360]), xlim([-5 365])
            %adf
            axis tight
            yL = get(gca,'ylim');
            if yL(1)>0
                yL(1) = 0;
            end
            yL(2) = yL(2)+.1;
            ylim(yL)

        end
    end
end


%%

trialdom_build = trialdom;
trialdom_predict = 1;
%revCorrPrediction(cellS.cellMat,trialdom_build,trialdom_predict)


% if ~isempty(idExamp)
%     for c = 1:length(DM.colordom)
%         figure
%         Nex = length(kern_examp{1});
%         for i = 1:Nex
%             
%             subplot(Nex,1,i)
%             tdom = 0:length(tcourse{c}{i})-1;
%             tdom = tdom*acqPeriod/1000;
%             plot(tdom,phi(tcourse{c}{i}))
%             %hold on
%            % plot(tdom,phi(Decon(phi(tcourse{c}{i}))),'k')         
%            
%             
%         end
%         xlabel('sec')
%     end
% end

%%
CH = GetTrialData([1 0],5);
im{1} = mean(CH{1}(:,:,2:end-1),3);
dum = im{1};

LP = fspecial('gaussian',size(im{1}),1); LP = LP/sum(LP(:));
HP = fspecial('gaussian',size(im{1}),5); HP = HP/sum(HP(:));
BP = LP-HP;
imH = ifft2(fft2(BP).*fft2(dum));
imH = fftshift(fftshift(imH,1),2);
W = 20;  
%%  This is for Eyal
%[xmicperpix ymicperpix] = getImResolution;
% 
% global maskSr
% 
% if ~isempty(idExamp)
%     for c = 1:length(DM.colordom)
%         figure
%         Nex = length(kern_examp{1});
%         for i = 1:Nex
%             %%%%%%%%%%%%%%%%           
%                                  
%             yrange = ((round(MK.CoM(idExamp(i),1)))-W/2:(round(MK.CoM(idExamp(i),1)))+W/2)-2;
%             xrange = (round(MK.CoM(idExamp(i),2)))-W/2:(round(MK.CoM(idExamp(i),2)))+W/2;            
%             subplot(Nex,4,i*4-3)
%             imagesc(imH(yrange,xrange)), colormap gray
%             hold on, 
%             se = strel('disk',1);
%             dum = imopen(maskS.bwCell{1}(yrange,xrange),se);
%             contour(dum,[.5 .5],'r')
%             
%             axis image         
%             axis off
%             
%             subplot(Nex,4,i*4-2)
%             imagesc(oridynamicsEx{c}{i}), %colormap jet
% %             try
% %             xl = str2num(cell2mat(get(gca,'XtickLabel')));
% %             catch
% %                'h' 
% %             end
% %             set(gca,'XtickLabel',round(DM.taudom(xl)))
% %             yl = str2num(cell2mat(get(gca,'YtickLabel')));
% %             set(gca,'YtickLabel',round(DM.oridom(yl)*10)/10,'TickDir','out'), axis square
%             xlabel('time'), ylabel('ori')
% 
%             
%             [idori idtau] = find(oridynamicsEx{c}{i} == max(oridynamicsEx{c}{i}(:)));
%             
%             subplot(Nex,4,i*4-1)
%             tcoursema = oridynamicsEx{c}{i}(idori,:)';
%             d90 = 90/(DM.oridom(2)-DM.oridom(1));
%             try
%                 tcourseOrth = oridynamicsEx{c}{i}(idori+d90,:)';
%             catch
%                 tcourseOrth = oridynamicsEx{c}{i}(idori-d90,:)';
%             end
%             
%             tsigma = tsigcoursema_examp{c}{i}/sqrt(3);
%             fill([DM.taudom fliplr(DM.taudom)],[tcoursema-tsigma; flipud(tcoursema+tsigma)]',[.0 .0 1])
%             hold on
%             plot(DM.taudom,tcoursema,'k'),xlabel('ms')
%             %ylim([min(tcoursema_examp{c}{i}-tsigcoursema_examp{c}{i})-.1 max(tcoursema_examp{c}{i}+tsigcoursema_examp{c}{i})+.1])
% 
%             tsigmi = tsigcoursemi_examp{c}{i}/sqrt(3);
%             fill([DM.taudom fliplr(DM.taudom)],[tcourseOrth-tsigmi; flipud(tcourseOrth+tsigmi)]',[1 .0 .0])
%             hold on
%             plot(DM.taudom,tcourseOrth,'k'),xlabel('ms')
%             xlim([DM.taudom(1) DM.taudom(end)])
%             legend('peak','','orthogonal','')
%             
%             [dum id] = max(var(oridynamicsEx{c}{i}));
%             oritcATpeak = oridynamicsEx{c}{i}(:,id);
%             [param ffit_ori varacc_ori sigma] = Gaussfit(DM.oridom,oritcATpeak',1);  param(2) = sigma;
%             
%             subplot(Nex,4,i*4)
%             plot(DM.oridom,oritcATpeak,'ok'),xlabel('ori'), title(['FWHM ~ ' num2str(round(sigma*.61*2))])
%             %hold on, plot(DM.oridomI,ffitoriEx{c}{i},'r')
%             set(gca,'Xtick',[0 90 180]), xlim([-5 180])
%             %ylim([min(ffitoriEx{c}{i})-.1 max(ffitoriEx{c}{i})+.1])
%             
% 
%         end
%     end
% end


function smoother = getSmoother(tausig,orisig,sfsig,taudom,oridom,sfdom)

if length(sfdom) == 1
    ksf = 1;
end
if length(oridom) == 1
    kori = 1;
end

ktau = exp(-(taudom-mean(taudom)).^2/(2*tausig^2)); %tausig is in ms

dom = log2(sfdom)-mean(log2(sfdom));
ksf = exp(-dom.^2/(2*sfsig^2)); %sfsig is in octaves

oridomdum = linspace(0,180,length(oridom)+1);
oridomdum = oridomdum(1:end-1);  %I use this one in case oridom wraps around
kori = exp(-(oridomdum-oridomdum(ceil(end/2))).^2/(2*orisig^2)); %orisig is in degrees

kdum = kori'*ksf;
smoother = zeros(length(kori),length(ksf),length(ktau));
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));


function tcori = intCirc(oridom,tcori,oridomI)

Ifac = length(oridomI)/length(oridom);

[dum idma] = max(tcori);
shifta = round(length(tcori)/2)-idma;
tcori = circshift(tcori(:)', [0 shifta]);

tcori = interp1(oridom,tcori,oridomI,'spline');

shiftb = round(-shifta*Ifac);
tcori = cir