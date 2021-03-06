function [kern kernblank countmat countmatblank] = Ggetrevcorrkernel(cellMat,tauN,trialdom,hh)

%Keep cellMat as an input because it filters it, aand I don't want to
%create a new variable of the same size

%3 takes the cell time courses as input, instead of all the images. 

global ACQinfo Analyzer cellS maskS G_RChandles

%%%%

blankNorm = get(G_RChandles.blankNorm,'value');

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
dtau = acqPeriod;
Ntau = round(tauN/acqPeriod);
taudom = 0:dtau:dtau*Ntau;


hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

blankProb = getparam('blankProb');

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});

sfdom_dum = domains.sfdom;
%Insert a place-holder for the blanks... the sequence will have index
%values that are one longer than the length of the spatial frequency
%domain, which are the blanks.
if blankProb > 0
    sfdom_dum = [sfdom_dum NaN];   
end


Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 
Tupdate = Tf*hper;


for s = 1:length(rseeds)
    
    eval(['colorS = rseed' num2str(s) '.colorseq;']);
    eval(['phaseS = rseed' num2str(s) '.phaseseq;']);
    eval(['sfS = rseed' num2str(s) '.sfseq;']);
    eval(['oriS = rseed' num2str(s) '.oriseq;']);    
    
    colorseq{s} = domains.colordom(colorS);   
    phaseseq{s} =  domains.phasedom(phaseS);  
    sfseq{s} =  sfdom_dum(sfS);     
    oriseq{s} =  domains.oridom(oriS);  
    
    %insert NaN for blanks
    if blankProb > 0
        idb = find(sfS == length(domains.sfdom)+1);
        
        colorseq{s}(idb) = NaN;
        phaseseq{s}(idb) =  NaN;
        sfseq{s}(idb) =  NaN;  %this is redundant
        oriseq{s}(idb) =  NaN;
    end
 
end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;


if ~isempty(hh)
    hh = ones(length(cellMat{1}(:,1)),1,'single')*hh(:)';
    
%     dom = (1:length(cellMat{1}(1,:)))';
%     H = [dom dom.^2 ones(length(dom),1)];
%     mat = inv(H'*H)*H';
    for T = 1:length(cellMat)       
        [c r] = getcondrep(T);
        
        %mu = mean(tcourse);
        %cellMat{T} = LFPfilt(cellMat{T},0,1000/acqPeriod,4,.05);
        
        %cellMat{T} = zscore(cellMat{T});
        
        %figure,plot(mean(cellMat{1}))
        
        %cellMat{T} = cellMat{T} - ones(length(cellMat{T}(:,1)),1)*cellMat{T}(1,:);
        
        
        cellMat{c}(:,:,r) = ifft(hh.*fft(squeeze(cellMat{c}(:,:,r)),[],2),[],2);

        
%         for k = 1:length(cellMat{T}(:,1))
%             dum = cellMat{T}(k,:)';
%             params = mat*dum;
%             dumhat = params(1)*H(:,1) + params(2)*H(:,2);
%             
%             cellMat{T}(k,:) = dum - dumhat;
% 
%         end            

    end
    
    
end


Ncell = length(celldom);
NT = getnotrials;
tcoursedum = 0;
figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    countmatblank{p} = zeros(1,length(taudom));
    kernblank{p} = zeros(1,length(taudom));


    for trialid = 1:length(trialdom)
        
        T = trialdom(trialid);
        [cond rep] = getcondrep(T);
        
%         dsync = diff(cellS.synctimes{cond,rep});
%         if any(abs(dsync-Tupdate/1000)>.100)
%             'Warning: syncs may be messed up'
%         end
       
        tcourse = squeeze(cellMat{cond}(p,:,rep));
        
%        tcourse = LFPfilt(tcourse,0,1000/acqPeriod,2,.05);
%                 fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
%                 fdom = fdom(1:end-1);
%                 figure,plot(fdom,abs(fft(tcoursedum-mean(tcoursedum))))
     
        
        %tcourse = zscore(tcourse);
        
        tdom = (0:length(tcourse)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;
        
        %idbase = find(tdom_pix<cellS.synctimes{cond,rep}(1)*1000 | tdom_pix>(cellS.synctimes{cond,rep}(end)*1000+500));
        %bLine = mean(tcourse(idbase));
        
        %tcourse = (tcourse-bLine)/bLine;
        
        seedno = Analyzer.loops.conds{cond}.val{1};

        oriseqdum = oriseq{seedno};
        sfseqdum = sfseq{seedno};       
        phaseseqdum = phaseseq{seedno};  
        colorseqdum = colorseq{seedno};         
        

        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)
                for phase = 1:length(phasedom)
                    for color = 1:length(colordom)
                        
                        id = find(oriseqdum == oridom(ori) & sfseqdum == sfdom(sf) & phaseseqdum == phasedom(phase)& colorseqdum == colordom(color));
                        
                        if ~isempty(id)
                            
                            stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times
                            %stimes = cellS.synctimes{cond,rep}(id)*1000;

                            for i = 1:length(stimes)
                                
                                idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
                                if ~isempty(idx1)
                                    idx1 = idx1(1);
                                    tpiece = idx1:idx1+length(taudom)-1;
                                    
                                    if tpiece(1)>0 & tpiece(end)<length(tcourse)
                                        tcoursePiece = squeeze(tcourse(tpiece))';
                                        kern{p}(ori,sf,phase,color,:) = squeeze(kern{p}(ori,sf,phase,color,:)) + tcoursePiece;
                                        countmat{p}(ori,sf,phase,color,:) = countmat{p}(ori,sf,phase,color,:) + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        id = find(isnan(oriseqdum));
        stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times
        %stimes = cellS.synctimes{cond,rep}(id)*1000;
        
        for i = 1:length(stimes)
            
            idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
            idx1 = idx1(1);
            tpiece = idx1:idx1+length(taudom)-1;
            
            if tpiece(1)>0 & tpiece(end)<length(tcourse)
                kernblank{p} = kernblank{p} + tcourse(tpiece);
                countmatblank{p} = countmatblank{p} + 1;
            end
        end
            
        
    end

    id = find(countmat{p} == 0);
    if ~isempty(id)
        'missing'
    end
 
    kern{p} = kern{p}./countmat{p};
    kernblank{p} = kernblank{p}./countmatblank{p};

    if blankNorm
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)
                for phase = 1:length(phasedom)
                    for color = 1:length(colordom)
                        
                        kern{p}(ori,sf,phase,color,:) = (squeeze(kern{p}(ori,sf,phase,color,:)) - kernblank{p}(:));
                        
                    end
                end
            end
        end
    end
    
    kernplot = kern{p}(:,1:3,:,:,:);
    kernplot = mean(kernplot,3);  %mean across phase
    kernplot = mean(kernplot,4); %mean across color
    
    %kernblank = kern{p}(end,:);
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    kernplot = squeeze(nanmean(kernplot,2));  %mean across spatial freqency
    blankmat = ones(length(kernplot(:,1)),1)*kernblank{p}(:)';

    %kernplot = (kernplot-blankmat);
    %kernplot = [kernplot; kernblank{p}(:)'];

    imagesc(real(kernplot))
    
    drawnow
end


