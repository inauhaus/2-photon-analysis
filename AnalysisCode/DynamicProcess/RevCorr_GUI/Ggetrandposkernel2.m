function [kern kernblank countmat countmatblank kernsig] = Ggetrandposkernel2(cellMat,trialdom,hh)

%2 gets the color info 

%Keep cellMat as an input because it filters it, and I don't want to
%create a new variable of the same size


global ACQinfo Analyzer cellS G_RChandles maskS

%%%%

tauN = str2num(get(G_RChandles.kernelLength,'string'));

%%%Get rid of the glia:

nID = getNeuronMask;  %get the index values for the neurons
%nID = 1;

masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);

%%%%

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)


%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);

hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

%expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
%load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});


logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt])


Tf = 1000/frate;  %Frame period in ms (frate obtained from log file)
Tupdate = Tf*hper;



%The appendage to 'rseed', e.g. 'rseed6', is the trial no., not the seed
%number.  I was making a mistake before.
Ntrials = getnotrials;
for t = 1:Ntrials  
    eval(['oriS = rseed' num2str(t) '.oriseq;']);
    eval(['xS = rseed' num2str(t) '.xseq;']);
    eval(['yS = rseed' num2str(t) '.yseq;']);
    eval(['bwS = rseed' num2str(t) '.bwseq;']);    
    eval(['colorS = rseed' num2str(t) '.colorseq;']);
    
    oriseq{t} =  domains.oridom(oriS);
    xseq{t} = domains.xdom(xS);   
    yseq{t} =  domains.ydom(yS);  
    bwseq{t} =  domains.bwdom(bwS);   
    colorseq{t} =  domains.colordom(colorS);
end


%used to do this
% for s = 1:length(rseeds)
% 
%     if exist(['rseed' num2str(s)])
%         eval(['bwS = rseed' num2str(s) '.bwseq;']);
%         eval(['xS = rseed' num2str(s) '.xseq;']);
%         eval(['yS = rseed' num2str(s) '.yseq;']);
%         eval(['oriS = rseed' num2str(s) '.oriseq;']);
%         eval(['colorS = rseed' num2str(s) '.colorseq;']);
% 
%         bwseq{s} = domains.bwdom(bwS);
%         xseq{s} =  domains.xdom(xS);
%         yseq{s} =  domains.ydom(yS);
%         oriseq{s} =  domains.oridom(oriS);
%         colorseq{s} =  domains.colordom(colorS);
%         
%         bwseq{s} = bwseq{s}(:);
%         xseq{s} =  xseq{s}(:);
%         yseq{s} =  yseq{s}(:);
%         oriseq{s} =  oriseq{s}(:);
%         colorseq{s} =  colorseq{s}(:);
%         
%     end
% 
% end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
xdom = domains.xdom;
ydom = domains.ydom;
bwdom = domains.bwdom;
colordom = domains.colordom;

Ncell = length(nID);
NT = getnotrials;
tcoursedum = 0;
figure
for p = 1:Ncell
    
    pID = nID(p);
    
    [idcelly idcellx] = find(masklabel == celldom(p));
       
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    
    countmat{p} = zeros(length(oridom),length(xdom),length(ydom),length(bwdom),length(colordom),length(taudom));
    kern{p} = zeros(length(oridom),length(xdom),length(ydom),length(bwdom),length(colordom),length(taudom));
    kernsq{p} = zeros(length(oridom),length(xdom),length(ydom),length(bwdom),length(colordom),length(taudom));
    countmatblank{p} = zeros(1,length(taudom));
    kernblank{p} = zeros(1,length(taudom));
    
    for trialid = 1:length(trialdom)
    %for trialid = 11:40

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);
        
%         dsync = diff(cellS.synctimes{cond,rep});
%         if any(abs(dsync-Tupdate/1000)>100)
%             'Warning: syncs may be messed up'
%         end                       

        tcourse = squeeze(cellMat{cond}(pID,:,rep));
        
        if p == 1
            %tcourse = squeeze(cellMat{cond}(1,:,rep)); %neuropil
            
            %tcourse = mean(squeeze(cellMat{cond}(:,:,rep))); %everything
            
            tcourse = mean(squeeze(cellMat{cond}(1:end,:,rep))); %cells
        end
        
        if ~isempty(hh)
            dsize = length(hh)-length(tcourse);
            if dsize>0
                tcourse = [tcourse ones(1,dsize)*tcourse(end)];
            elseif dsize<0
                tcourse = tcourse(1:end+dsize);
            end
        end
        
        
        tcourse = processTcourse(tcourse,hh,1,acqPeriod);
        
        %        tcourse = LFPfilt(tcourse,0,1000/acqPeriod,4,.05);
        %                 fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
        %                 fdom = fdom(1:end-1);
        %                 figure,plot(fdom,abs(fft(tcoursedum-mean(tcoursedum))))        
        
        %tcourse = zscore(tcourse);
        
        tdom = (0:length(tcourse)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;
        
        %idbase = find(tdom_pix<synctimes{cond,rep}(1)*1000 | tdom_pix>(synctimes{cond,rep}(end)*1000+500));
        %bLine = mean(tcourse(idbase));
        
        %tcourse = (tcourse-bLine)/bLine;        
        
        
        %This is not correct all the time... only if trial no matches seed
        %no.       
%         seedno = Analyzer.loops.conds{cond}.val{1};        
%         oriseqdum = oriseq{seedno};
%         xseqdum = xseq{seedno};
%         yseqdum = yseq{seedno};
%         bwseqdum = bwseq{seedno};
%         colorseqdum = colorseq{seedno};

        oriseqdum = oriseq{T}(:);
        xseqdum = xseq{T}(:);       
        yseqdum = yseq{T}(:);  
        bwseqdum = bwseq{T}(:); 
        colorseqdum = colorseq{T}(:);
        
        for ori = 1:length(oridom)
            for x = 1:length(xdom)
                for y = 1:length(ydom)
                    for bw = 1:length(bwdom)
                        for color = 1:length(colordom)

                            
                            if size(oriseqdum(:),2) == 1
                                id = find(oriseqdum == oridom(ori) & xseqdum == xdom(x) & yseqdum == ydom(y) & bwseqdum == bwdom(bw) & colorseqdum == colordom(color));                                                               
                            else
                                id = [];
                                for drop = 1:size(size(oriseqdum,2))
                                    iddum = find(oriseqdum(:,drop) == oridom(ori) & xseqdum(:,drop) == xdom(x) & yseqdum(:,drop) == ydom(y) & bwseqdum(:,drop) == bwdom(bw) & colorseqdum(:,drop) == colordom(color));
                                    id = [id(:); iddum(:)];
                                end
                            end
                            
                            %stimes = (id-1)*Tupdate; %Stimulus times
                            %stimes = cellS.synctimes{cond,rep}(id)*1000;
                            stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                            for i = 1:length(stimes)

                                [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1))));
                                %idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
                                idx1 = idx1(1);
                                tpiece = idx1:idx1+length(taudom)-1;

                                if tpiece(1)>0 & tpiece(end)<length(tcourse)
                                    kern{p}(ori,x,y,bw,color,:) = squeeze(kern{p}(ori,x,y,bw,color,:)) + squeeze(tcourse(tpiece))'; %E(x)
                                    kernsq{p}(ori,x,y,bw,color,:) = squeeze(kernsq{p}(ori,x,y,bw,color,:)) + squeeze(tcourse(tpiece))'.^2; %E(x^2)
                                    countmat{p}(ori,x,y,bw,color,:) = countmat{p}(ori,x,y,bw,color,:) + 1;
                                end

                            end
                        end
                    end
                end
            end
        end
        
        
    end
    
    %Reminder:  randpos does not have any blanks
    
    if length(ydom) > 1        
        kern{p} = sum(kern{p},3);
        kernsq{p} = sum(kernsq{p},3);
        countmat{p} = sum(countmat{p},3);        
    end
    
    kern{p} = kern{p}./countmat{p};  
    kernsq{p} = kernsq{p}./countmat{p};    
    kernsig{p} = sqrt(kernsq{p} - kern{p}.^2);  %E((x-u)^2) = E(x^2) - (E(x))^2    

    
    kern{p} = reshape(kern{p},[length(oridom) length(xdom) length(bwdom) length(colordom) length(taudom)]);  %get rid of 'y' dimension
    kernsig{p} = reshape(kernsig{p},[length(oridom) length(xdom) length(bwdom) length(colordom) length(taudom)]);  %get rid of 'y' dimension
    kerncount{p} = reshape(countmat{p},[length(oridom) length(xdom) length(bwdom) length(colordom) length(taudom)]);  %get rid of 'y' dimension
    
    
    cellS.kernAll = kern;
    cellS.kernSigAll = kernsig;
    cellS.kernCount = kerncount;
    
    %Downsample
    kernplot = kern{p};
       
    countplot = countmat{p};
    idmissing = find(countplot(:) == 0);

    %if ~isempty(idmissing)
        if(1)
        
        kernplot(idmissing) = 0;
        
        Neven = size(kernplot,2) - rem(size(kernplot,2),2);
        dumA = kernplot(:,1:2:Neven,:,:,:,:); dumB = kernplot(:,2:2:Neven,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);
        
        kernplot = (dumA + dumB)/2; %oridomain
        xdomplot = (xdom(1:2:Neven) + xdom(2:2:Neven))/2;
        
        
        countplot = countplot(:,1:2:Neven,:,:,:,:) + countplot(:,2:2:Neven,:,:,:,:);
        
        Neven = size(kernplot,1) - rem(size(kernplot,1),2);
        dumA = kernplot(1:2:Neven,:,:,:,:,:); dumB = kernplot(2:2:Neven,:,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);

        kernplot = (dumA + dumB)/2; %spatial domain
        oridomplot = (oridom(1:2:Neven) + oridom(2:2:Neven))/2;
        
        countplot = countplot(1:2:Neven,:,:,:,:,:) + countplot(2:2:Neven,:,:,:,:,:);
        

        'Downsampling because not enough presentations. Only for plotting here'

    end
    
    vartime = mean(mean(mean(mean(kernplot.^2,1),2),3),4);
    vartime = squeeze(sqrt(vartime));
    
    
    
    idt = find(taudom>50 & taudom<400);
    
    kernplot = mean(kernplot(:,:,:,:,idt),5);  %mean across time
    kernplot = mean(kernplot(:,:,:,:,:),4);  %mean across color
    
    kernplot = squeeze(kernplot); 

    kernplot_b = squeeze(kernplot(:,:,1,:));
    kernplot_w = squeeze(kernplot(:,:,2,:));
    
    %kernblank = kern{p}(end,:);

    
    kernsmooth = zeros(size(kernplot_b));
    kernsmooth(1:3,1:3) = [.2 1 .2]'*[.2 1 .2];
    kernplot_b = ifft2(fft2(kernplot_b).*abs(fft2(kernsmooth)));
    kernplot_w = ifft2(fft2(kernplot_w).*abs(fft2(kernsmooth)));
     kernplot = kernplot_b'+kernplot_w';
    
%     dum = mean(kernplot_b,3);
%     [idy idx] = find(dum == max(dum(:)))
%     ttrace = squeeze(kernplot_b(idy,idx,:));
%     
     subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
%     kernplot =  kernplot_w'+kernplot_b';
%     kernplot = kernplot-prctile(kernplot(:),10);
    
    %kernplot(find(kernplot<0)) = 0;
    RF = iradon(kernplot,oridomplot);
    
    %imagesc(iradon([zeros(5,18); kernplot'; zeros(5,18)],oridom))
    imagesc(kernplot)
    %imagesc(RF)
    %plot(ttrace)
    %ylim([-.5 1])
    %plot(taudom,vartime)
    
    %plot([mean(kernplot_b(5:7,:))' mean(kernplot_w(5:7,:))'])
    
    drawnow
end


