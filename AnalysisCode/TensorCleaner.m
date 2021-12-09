function CH = TensorCleaner(CH,slowMoFlag,fastXcorrFlag,opticFlowFlag,frameRemovalFlag,SVDremovalFlag,DSflag,Fwin,varargin)        

global cellS maskS
    
Atemp = maskS.anatomyTemplate;


if ~isempty(varargin)
    
    paramVec = varargin{1};
    
    opticFlowXsig = paramVec(1); %smooth in X before computing optic flow
    rthresh_floor = paramVec(2); %absolute floor; correlation coefficient with the template.
    deletionThresh_Z = paramVec(3); %Threshold below the median in SD
    N_PC = paramVec(4); %number of principal components to remove
else %defaults
    opticFlowXsig = 2; %smooth in X before computing optic flow
    rthresh_floor = .3; %absolute floor; correlation coefficient with the template.
    deletionThresh_Z = 1.5; %Threshold below the median in SD
    N_PC = 1; %number of principal components to remove
end

   

if slowMoFlag %Align the entire block by a single dx dy shift. 
    
    imdum = median(CH,3); %align to the median   
    
    [mbestSlow nbestSlow] = getShiftVals(imdum.^2,Atemp.^2,[0 0]);  %squaring seems to really help sometimes
    
    CH = circshift(CH,[-round(mbestSlow) -round(nbestSlow) 0]);
    
end




if fastXcorrFlag
    

    %% Coarse x/y aligmnment using the motion info in cellS. Uses template match.
    
    
    if Fwin(2) > length(cellS.motionInfo.mbestFrame)  %This handles the last trial from CondMaskdata2
        mbest = cellS.motionInfo.mbestFrame(Fwin(1):end);
        nbest = cellS.motionInfo.nbestFrame(Fwin(1):end);
        app = size(CH,3)-length(mbest); %append
        mbest = [mbest mbest(end)*ones(1,app)];
        nbest = [nbest nbest(end)*ones(1,app)];
    else
        mbest = cellS.motionInfo.mbestFrame(Fwin(1):Fwin(2));
        nbest = cellS.motionInfo.nbestFrame(Fwin(1):Fwin(2));
    end

    
    mbest = medfilt1(mbest,3); %remove jitter
    nbest = medfilt1(nbest,3);
    CH = align2Template(CH,mbest,nbest);
    
end

if opticFlowFlag
    
    %% Shift correction using optic flow. Refinement.
    
    %sigx = [3 1];
    
    [CH dx dy] = TensorOpticFlowRigidCorrection(CH, opticFlowXsig/(DSflag+1));
    
    
    %% Get motion information and subtract it [using the unsmoothed version]
    
    %[CH dx dy] = TensorSubtractFlow(CH);
    
    
end


%% Bin each image in space and time
 

if DSflag  
        
    CH = Tensor_2x2bin(CH); %2x2 binning
    
    % Smooth in time, then downsample
%     Lsig = 1;
%     Hsig = inf;
%     CH = TensorFilter_time(CH,Lsig,Hsig); %smooth in time;
    %CH = CH(:,:,1:D:end);
    
%     rema = rem(size(CH,3),2);
%     CH = CH(:,:,1:2:(end-rema)) + CH(:,:,2:2:(end-rema));
    
end
CHdim = size(CH);


%%  Subtract residual motion

% if fastMoFlag
% 
%     [CH dx dy] = TensorSubtractFlow(CH);
% 
% end

%%%%%%%


%%
idGoodFrames = ones(1,CHdim(3)); %default
if frameRemovalFlag  %% NaN bad frames (before SVD removal below)
    
    rTempAll = cellS.motionInfo.rTempFrames;
    
    if Fwin(2) > length(rTempAll)  %This handles the last trial from CondMaskdata2
        rTemp = rTempAll(Fwin(1):end);
        app = size(CH,3)-length(rTemp); %append
        rTemp = [rTemp rTemp(end)*ones(1,app)];
    else
        rTemp = rTempAll(Fwin(1):Fwin(2));
    end
    
           
    %Get threshold based on stats of the entire experiment, not just the
    %trial.
    [rthresh] = getrThresh(rTempAll,deletionThresh_Z);
    idGoodFrames(find(rTemp<rthresh)) = 0;  %Find bad frames within this window of the experiment
    
    CH = reshape(CH,[CHdim(1)*CHdim(2) CHdim(3)])'; %vectorize
        
    CH(find(1-idGoodFrames),:) = NaN;
    
    %[CH rTemp] = TensorNaNBadFrames(CH,maskS.anatomyTemplate,rthresh,rthresh_floor);
    
    CH_NaNdump = CH;
    CH_NaNdump(find(1-idGoodFrames),:) = [];
    
    if sum(idGoodFrames)<0.5*CHdim(3)  %abandon this trial if it totally sucks
        CH = CH*NaN;
        idGoodFrames = [];    
    end
    
else 
    CH_NaNdump = CH;
    
end

if SVDremovalFlag & sum(idGoodFrames)
    
    if length(size(CH)) == 3
        CH = reshape(CH,[CHdim(1)*CHdim(2) CHdim(3)])'; %vectorize
    end
    
    if length(size(CH_NaNdump)) == 3
        CH_NaNdump = reshape(CH_NaNdump,[CHdim(1)*CHdim(2) CHdim(3)])'; %vectorize
    end
    
    CH_NaNdump = TensorRemoveSingularValues(CH_NaNdump,N_PC);
    
    CH(find(idGoodFrames),:) = CH_NaNdump; %Replace with the PCA subtracted matrix. The others are NaN;
    
end

%Unvectorize
if length(size(CH)) == 2
    CH = reshape(CH',[CHdim(1) CHdim(2) size(CH,1)]);
end



% %%
% mi = prctile(CH(:),10)
% ma = prctile(CH(:),99.9)
% %
% 
% for i =1 : size(CH,3)
%     
%     figure(20)
%    imagesc((CH(:,:,i)),[mi ma]), colorbar
%    pause(.001)
%    clf
% end


