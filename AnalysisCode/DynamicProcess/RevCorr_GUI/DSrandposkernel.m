function [kernAll kernCount] = DSrandposkernel(kernAll,kernCount,DSx,DSo,DSBW)

countA = 0; countB = 0;
for p = 1:length(kernAll)
    
        id = find(isnan(kernAll{p}));
        kernAll{p}(id) = 0; %Turn NaN spots to zeros to make stuff below work!
        
    if DSBW
        
        %First combine black and white kernels:
        kernB = kernAll{p}(:,:,1,:,:,:);
        countB = kernCount{p}(:,:,1,:,:,:);

        kernA = kernAll{p}(:,:,2,:,:,:);
        countA = kernCount{p}(:,:,2,:,:,:);
        
        kernAll{p} = kernA.*countA + kernB.*countB; %"un-normalize"
        kernAll{p} = kernAll{p}./(countA + countB);  %downsample and renormalize   
        
        %kernAll{p} = kernA; 
 
        
        id = find((countA + countB) == 0);
        kernAll{p}(id) = 0;
        kernCount{p} = countA + countB;
    end
    
    if DSx>1

        id = find(isnan(kernAll{p}));
        kernAll{p}(id) = 0;
        %

        re = rem(size(kernAll{p},2),2);
        kernB = kernAll{p}(:,2:2:end,:,:,:,:);
        countB = kernCount{p}(:,2:2:end,:,:,:,:);

        kernA = kernAll{p}(:,1:2:end-re,:,:,:,:);
        countA = kernCount{p}(:,1:2:end-re,:,:,:,:);
        
%         kernA = kernAll{p};
%         countA = kernCount{p};
%         
%         kernB = circshift(kernAll{p},[0 1 0 0 0 0]);
%         countB = circshift(kernCount{p},[0 1 0 0 0 0]);

        kernAll{p} = kernA.*countA + kernB.*countB; %"un-normalize"
        kernAll{p} = kernAll{p}./(countA + countB);  %downsample and renormalize
        kernCount{p} = countA + countB;

%         kernAll{p} = kernAll{p}(:,1:2:end,:,:,:,:);
%         kernCount{p} = kernCount{p}(:,1:2:end,:,:,:,:);

    
        id = find((countA + countB) == 0);
        kernAll{p}(id) = 0;

    end

    if DSo

        id = find(isnan(kernAll{p}));
        kernAll{p}(id) = 0;
% 
        kernA = kernAll{p}(1:2:end,:,:,:,:,:);
        countA = kernCount{p}(1:2:end,:,:,:,:,:);

        kernB = kernAll{p}(2:2:end,:,:,:,:,:);
        countB = kernCount{p}(2:2:end,:,:,:,:,:);
        
%         kernA = kernAll{p};
%         countA = kernCount{p};
%         
%         kernB = circshift(kernAll{p},[1 0 0 0 0 0]);
%         countB = circshift(kernCount{p},[1 0 0 0 0 0]);

        kernAll{p} = kernA.*countA + kernB.*countB; %"un-normalize"
        kernAll{p} = kernAll{p}./(countA + countB);  %downsample and renormalize
        kernCount{p} = countA + countB;

%         kernAll{p} = kernAll{p}(2:2:end,:,:,:,:,:);
%         kernCount{p} = kernCount{p}(2:2:end,:,:,:,:,:);
        
%         id = find((countA + countB) == 0);
%         kernAll{p}(id) = 0;


    end

end
