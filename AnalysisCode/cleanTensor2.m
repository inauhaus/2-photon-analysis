function ToI = cleanTensor2(T,Fwin)

%Will return an array that is downsampled 2x

global maskS cellS G_handles

%%x/y aligmnment

fastMo = get(G_handles.fastMotionFlag,'Value');  %lateral movement correction

frames = Fwin(1):Fwin(end);

if fastMo
    
    if isfield(cellS.motionInfo,'mbest')
        
        for q = 1:size(T,3) %loop each frame
            mbest = cellS.motionInfo.mbestFrame(frames(q));
            nbest = cellS.motionInfo.nbestFrame(frames(q));
            
            %T(:,:,q) = circshift(T(:,:,q),[round(-mbest) round(-nbest) 0]);            
            T(:,:,q) = circshift_continous2(T(:,:,q),-nbest,-mbest);
            
            %                     if cellS.motionInfo.rTemp{t}(q) < rthresh
            %                         CH(:,:,q) = NaN;
            %                     end
        end
    else
        warning('Not applying motion correction.  Need to first "Get Motion Information".')
    end
end

%%

Atemp = maskS.anatomyTemplate;  %This could be from a different experiment if desired when setting directory
A = Atemp(1:2:end,1:2:end,:);
B = Atemp(1:2:end,2:2:end,:);
C = Atemp(2:2:end,1:2:end,:);
D = Atemp(2:2:end,2:2:end,:);
Atemp = A+B+C+D;
Atemp = Atemp(:);
Atemp = Atemp/norm(Atemp);

A = T(1:2:end,1:2:end,:);
B = T(1:2:end,2:2:end,:);
C = T(2:2:end,1:2:end,:);
D = T(2:2:end,2:2:end,:);
To = A+B+C+D;
dim = size(To);
To = reshape(To,[dim(1)*dim(2) dim(3)])';
dimV = size(To);


% h = 1./fft(mean(imR,2));
% h = h*ones(1,size(imR,2));
% imR = ifft(abs(fft(h)).*fft(imR));

%Smooth each pixel in time
tsig = 1;
if tsig>0
    h = fspecial('gaussian',[size(To,1) 1],tsig);
    h = h*ones(1,size(To,2));
    To = ifft(abs(fft(h)).*fft(To));
end
ToI = To(1:2:end,:);

muTo = mean(ToI,2)*ones(1,dimV(2));
sigTo = std(ToI,[],2)*ones(1,dimV(2));
ToI = (ToI-muTo)./sigTo;

Atemp = (Atemp-mean(Atemp(:)))/std(Atemp(:));
rTemp = (ToI*Atemp)/size(ToI,2); %Correlation coef with the template, at each frame

%ToI = subtractMotion(ToI,rTemp);





%% Remove bad frames

figure, plot(rTemp)

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

%Remove first singular value


ToI = ToI - ones(size(ToI,1),1)*median(ToI); %Center the "cloud"
ToI = ToI./(ones(size(ToI,1),1)*std(ToI)); %Center the "cloud"


[U,S,V] = svd(ToI','econ');
C = U(:,1)*V(:,1)'*S(1,1) ;
C = C + U(:,2)*V(:,2)'*S(2,2) + U(:,3)*V(:,3)'*S(3,3) ;
ToI = ToI-C';




%%

ToI = reshape(ToI',[dim(1) dim(2) size(ToI,1)]);
%%

% mi = prctile(ToI(:),10);
% ma = prctile(ToI(:),99.9);
% 
% %ToI = reshape(ToI',[dim(1) dim(2) size(ToI,1)]);
% for i =1 : size(ToI,3)
%     
%     figure(20)
%    imagesc((ToI(:,:,i)),[mi ma]), colorbar
%    pause(.2)
% end


% dim = size(imXdum);
% dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
% imXdum = interp1(1:dim(1),imXdum,linspace(1,dim(1),dimI(1)));
% imXdum = interp1(1:dim(2),imXdum',linspace(1,dim(2),dimI(2)))';





