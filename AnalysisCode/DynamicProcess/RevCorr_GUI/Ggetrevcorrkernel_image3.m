function [kernIm kernblank countmat countmatblank] = Ggetrevcorrkernel_image3(trialdom,hh)

%Ian Nauhaus

%2 incorporates new motion correction stuff
%3 incorporates the newest motion correction stuff, including SVD and frame
%deletion using the function TensorCleaner

%Keep cellMat as an input because it filters it, and I don't want to
%create a new variable of the same size

%3 takes the cell time courses as input, instead of all the images. 

global ACQinfo Analyzer G_RChandles G_handles domains cellS

%%%%

blankNorm = get(G_RChandles.blankNorm,'value');

%%%%

%%%Get the flags for image movement correction%%%%

slowMoFlag = get(G_handles.slowMotionFlag,'Value');  %lateral movement correction
fastMoFlag = get(G_handles.fastMotionFlag,'Value');  %lateral movement correction

frameRemovalFlag = get(G_handles.frameRemoval,'Value');  %lateral movement correction
deletionThresh_Z = str2num(get(G_handles.deletionThreshold,'string'));

SVDremovalFlag = get(G_handles.PCremoval,'Value');  %lateral movement correction
N_PC = str2num(get(G_handles.nPC,'string'));

opticFlowXsig = 2; %smooth in X before computing optic flow
rthresh_floor = .3; %absolute floor; correlation coefficient with the template.
%deletionThresh_Z = 1.5; %Threshold below the median in SD
%N_PC = 1; %number of principal components to remove

paramVec = [opticFlowXsig rthresh_floor deletionThresh_Z N_PC];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)


hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

blankProb = getparam('blankProb');

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

[domains seqs] = getSeqInfo(trialdom);

oridom = domains{4}.oridom;
sfdom = domains{4}.sfdom;
phasedom = domains{4}.phasedom;
colordom = domains{4}.colordom;


Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 
Tupdate = Tf*hper;

NT = getnotrials;

tau_xy = ACQinfo.linesPerFrame/2*ACQinfo.msPerLine;

countmat = cell(length(oridom),length(sfdom),length(phasedom),length(colordom));
kernIm = cell(length(oridom),length(sfdom),length(phasedom),length(colordom));
% countmatblank = zeros(1,length(taudom));
% kernblank = zeros(1,length(taudom));

trialstartID = find(ACQinfo.SBInfo.event_id == 3); %rising edge ID
trialendID = find(ACQinfo.SBInfo.event_id == 2); %falling edgeID
FramesPerTrial = min(ACQinfo.SBInfo.frame(trialendID) - ACQinfo.SBInfo.frame(trialstartID));

NT = getnotrials;


alignCh = 1;
chvec = [0 0 0 0];
chvec(alignCh) = 1;
temp = GetTrialData(chvec,1);
Tdim = size(temp{1});


if ~isempty(hh)
    hfilt = ones(Tdim(1),Tdim(2),length(hh));
    for m = 1:length(hh)
        hfilt(:,:,m) = hh(m)*hfilt(:,:,m);
    end
end


for trialid = 1:length(trialdom)

    T = trialdom(trialid);

T
    [CHraw f1 fN] = GetTrialData([1 0],trialdom(trialid));

    %%
    
    %append or truncate data tensor to match length of shortest trial.
    CH = CHraw{1}(:,:,1:FramesPerTrial);
    fN = f1+FramesPerTrial-1;
    
    
    %% This takes a while
    CHdim = size(CH);
    
    CH = TensorCleaner(CH,slowMoFlag,fastMoFlag,frameRemovalFlag,SVDremovalFlag,[f1 fN],paramVec);
    
    idGoodFrames = ones(1,CHdim(3));
    idGoodFrames(find(isnan(CH(1,1,:)))) = 0;
    
    CH(:,:,find(1-idGoodFrames)) = 0; %Change them from NaN to zero
    
    %%
    mi = prctile(CH(:),1);
    ma = prctile(CH(:),99.8);
    figure(20)
    for i = 1:size(CH,3)
        cla
        imagesc(CH(:,:,i),[mi ma])
        pause(.1)
        drawnow
        
    end
    
%%

    
    T = trialdom(trialid);
    [cond rep] = getcondrep(T);

    %         dsync = diff(cellS.synctimes{cond,rep});
    %         if any(abs(dsync-Tupdate/1000)>.100)
    %             'Warning: syncs may be messed up'
    %         end

    N = length(CH(1,1,:));
    tdom = (0:N-1)*acqPeriod;
    %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
    tdom_pix = tdom + tau_xy;

    seedno = Analyzer.loops.conds{cond}.val{1};

    for ori = 1:length(oridom)
        for sf = 1:length(sfdom)
            for phase = 1:length(phasedom)
                for color = 1:length(colordom)

                    id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq  == sfdom(sf) & seqs{T}.phaseseq  == phasedom(phase)& seqs{T}.colorseq  == colordom(color));

                    if ~isempty(id)
                        
                        if isempty(kernIm{ori,sf,phase,color})
                            kernIm{ori,sf,phase,color} = 0;
                            countmat{ori,sf,phase,color} = 0;
                        end
                        stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times
                        %stimes = cellS.synctimes{cond,rep}(id)*1000;

                        for i = 1:length(stimes)
                            [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                            if ~isempty(idx1)
                                idx1 = idx1(1);
                                tpiece = idx1:idx1+length(taudom)-1;

                                if tpiece(1)>0 & tpiece(end)<N
                                    kernIm{ori,sf,phase,color} = kernIm{ori,sf,phase,color} + CH(:,:,tpiece) ;
                                    %kernIm{ori,sf,phase,color} = kernIm{ori,sf,phase,color} - CH(:,:,1);
                                    countmat{ori,sf,phase,color} = countmat{ori,sf,phase,color} + idGoodFrames(tpiece);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    id = find(isnan(seqs{T}.oriseq));
    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times

    %stimes = cellS.synctimes{cond,rep}(id)*1000;

%     for i = 1:length(stimes)
% 
%         idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
%         idx1 = idx1(1);
%         tpiece = idx1:idx1+length(taudom)-1;
% 
%         if tpiece(1)>0 & tpiece(end)<length(tcourse)
%             kernblank{p} = kernblank{p} + tcourse(tpiece);
%             countmatblank{p} = countmatblank{p} + 1;
%         end
%     end


end


for ori = 1:length(oridom)
    for sf = 1:length(sfdom)
        for phase = 1:length(phasedom)
            for color = 1:length(colordom)
                for t = 1:length(taudom)
                    kernIm{ori,sf,phase,color}(:,:,t) = kernIm{ori,sf,phase,color}(:,:,t)/countmat{ori,sf,phase,color}(t);
                end
                % kernblank = kernblank./countmatblank;
            end
        end
    end
end
kernblank = 0;
% kernblank{p} = kernblank{p}./countmatblank{p};
% 
% if blankNorm
%     for ori = 1:length(oridom)
%         for sf = 1:length(sfdom)
%             for phase = 1:length(phasedom)
%                 for color = 1:length(colordom)
% 
%                     kernIm{p}(ori,sf,phase,color,:) = (squeeze(kernIm{p}(ori,sf,phase,color,:)) - kernblank{p}(:));
% 
%                 end
%             end
%         end
%     end
% end

    function y = subtractTemporalPoly(x,order)
        
        %subtracts polynomial fit from each pixel time course.  x is
        %assumed to be [x y T] dimensional
        
        xdim = size(x);
        
        dum = shiftdim(x,2); %make time the first dimension
        dum = reshape(dum,[xdim(3) xdim(1)*xdim(2)]);  %make it a 2D matrix, with time as first dim
        
        dom = (0:size(dum,1)-1)';
        
        H = ones(length(dom),order+1);  %last column for DC
        for i = 1:order
            H(:,i) = dom.^i;
        end
        
        p = inv(H'*H)*H'*dum;
        xfit = H*p;
        
        xfit = reshape(xfit,[xdim(3) xdim(1) xdim(2)]);
        xfit = shiftdim(xfit,1);
        y = x-xfit;

