function ToI = cleanTensor4(T,Fwin,D)

%If D = 2, will return an array that is downsampled 2x. Happens just before
%SVD.

global maskS cellS G_handles


%% Coarse x/y aligmnment using the motion info in cellS. Uses template match.

frames = Fwin(1):Fwin(end);


if get(G_handles.fastXcorr,'value')
    
    if isfield(cellS,'motionInfo')
        mbest = cellS.motionInfo.mbestFrame(frames);
        nbest = cellS.motionInfo.nbestFrame(frames);
        T = align2Template(T,mbest,nbest);
    else
        a = questdlg('There is no motion information.  Proceed to get the mask w/o motion correction?');
        if strcmp('No',a)
            asdf
        end
    end
    
end

%% Optic flow based correction methods.

if get(G_handles.OpticFlowCorrection,'value')
    
    %sigx = [3 1];
    sigx = 2;
     
    %Compute x/y shift and apply it
    [To dx dy] = TensorOpticFlowRigidCorrection(T, sigx);
    
    % Get motion information and subtract it [using the unsmoothed version]    
    [To dx dy] = TensorSubtractFlow(To);
    
else
    
    To = T;
    
end


%% Bin each image in space and time
 
%Binning is automatic.  It does not use input from GUI

Atemp = maskS.anatomyTemplate;  %This could be from a different experiment if desired when setting directory

if D > 1  
        
    Atemp = Tensor_2x2bin(Atemp);  %2x2 binning
    Atemp = Atemp(:);
    
    To = Tensor_2x2bin(To); %2x2 binning
    
    % Smooth in time, then downsample
    Lsig = D-1;
    Hsig = inf;
    To = TensorFilter_time(To,Lsig,Hsig); %smooth in time;
    ToI = To;
    ToI = ToI(:,:,1:D:end);
    
else
    ToI = To;
    Atemp = Atemp(:);
    
end
dim = size(ToI);
ToI = reshape(ToI,[dim(1)*dim(2) dim(3)])'; %vectorize


%% Get bad frames

if get(G_handles.frameRemoval,'value')
    
    %Zthresh = 1; %number of SDs below median for removal
    
    Zthresh = str2num(get(G_handles.deletionThreshold,'string'));
    rFloor = 0.4;
    [ToI rTemp] = TensorNaNBadFrames(ToI,Atemp,Zthresh,rFloor);
    
    
end

%% SVD removal
if get(G_handles.PCremoval,'value')
    
    idNaN = find(isnan(ToI(:,1)));
    ToNaNdump = ToI;
    ToNaNdump(idNaN,:) = [];
    
    N = str2num(get(G_handles.nPC,'string')); %number of SVs to remove
    ToNaNdump = TensorRemoveSingularValues(ToNaNdump,N);
    
    idnum = find(~isnan(ToI(:,1)));
    ToI(idnum,:) = ToNaNdump;
    
end


%idgoodTimes = find(rTemp>thresh);

%ToI = ToI(idgoodTimes,:);

%%

ToI = reshape(ToI',[dim(1) dim(2) size(ToI,1)]);
%%

%mi = prctile(ToI(:),10);
%ma = prctile(ToI(:),99.9);
%%

% for i =1 : size(ToI,3)
%     
%     figure(20)
%    imagesc((ToI(:,:,i)),[mi ma]), colorbar
%    pause(.001)
%    clf
% end


% dim = size(imXdum);
% dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
% imXdum = interp1(1:dim(1),imXdum,linspace(1,dim(1),dimI(1)));
% imXdum = interp1(1:dim(2),imXdum',linspace(1,dim(2),dimI(2)))';





