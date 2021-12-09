function ToI = cleanTensor3(T,Fwin)

%Will return an array that is downsampled 2x

global maskS cellS G_handles


%% Coarse x/y aligmnment using the motion info in cellS. Uses template match.

frames = Fwin(1):Fwin(end);
mbest = cellS.motionInfo.mbestFrame(frames);
nbest = cellS.motionInfo.nbestFrame(frames);
T = align2Template(T,mbest,nbest); 

%% Shift correction using optic flow. Refinement.

%sigx = [3 1];
sigx = 2;
To = TensorOpticFlowRigidCorrection(T, sigx);


%% Get motion information and subtract it [using the unsmoothed version]

To = TensorSubtractFlow(To);

%% Bin each image

Atemp = maskS.anatomyTemplate;  %This could be from a different experiment if desired when setting directory

Atemp = Tensor_2x2bin(Atemp);  %2x2 binning
Atemp = Atemp(:);

To = Tensor_2x2bin(To); %2x2 binning

%% Smooth in time, then downsample

%Smooth each pixel in time, and downsample (for speed)
Lsig = 1;
Hsig = inf;
To = TensorFilter_time(To,Lsig,Hsig); %smooth in time;
ToI = To;
ToI = ToI(:,:,1:2:end);

%% Get bad frames

dim = size(ToI);
ToI = reshape(ToI,[dim(1)*dim(2) dim(3)])'; %vectorize

ToZ = TensorZscore_space(ToI);

Atemp = (Atemp-mean(Atemp(:)))/std(Atemp(:));
rTemp = (ToZ*Atemp)/size(ToZ,2); %Correlation coef with the template, at each frame


%figure, plot(rTemp)

thresh = prctile(rTemp,50) - std(rTemp)

%thresh = 0.6;
idgoodTimes = find(rTemp>thresh);

ToI = ToI(idgoodTimes,:);

%%

% Subtract polynomial fit from each pixel (i.e. HP filter)
% H = (1:size(ToI,1))';
% %H = [H H.^2 H.^3 ones(length(H(:,1)),1)];
% H = [H ones(length(H(:,1)),1)];
% slps = inv(H'*H)*H'*ToI;
% imRfit = H*slps;
% ToI = ToI-imRfit;



%%



N = 3; %number of SVs to remove
ToI = TensorRemoveSingularValues(ToI,N);


%%

ToI = reshape(ToI',[dim(1) dim(2) size(ToI,1)]);
%%

mi = prctile(ToI(:),10);
ma = prctile(ToI(:),99.9);
%%

for i =1 : size(ToI,3)
    
    figure(20)
   imagesc((ToI(:,:,i)),[mi ma]), colorbar
   pause(.001)
   clf
end


% dim = size(imXdum);
% dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
% imXdum = interp1(1:dim(1),imXdum,linspace(1,dim(1),dimI(1)));
% imXdum = interp1(1:dim(2),imXdum',linspace(1,dim(2),dimI(2)))';





