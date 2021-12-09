function [To dx dy] = TensorOpticFlowRigidCorrection(T, sigx)

To = T; %'To' gets iteratively corrected


for q = 1:length(sigx)
    
%     hh = fspecial('gaussian',size(T(:,:,1)),sigx(q));
%     hh = abs(fft2(hh));
%     TSmooth = zeros(size(T));
%     for i = 1:size(To,3)
%         im = To(:,:,i);
%         TSmooth(:,:,i) = ifft2(fft2(im).*hh);
%     end
%     
    TSmooth = TensorFilter_space(To,sigx(q),inf); %filter in space
    
    %imgs_uncorrected = smoothTensorTime(imgs_uncorrected,tsig); %smooth in time
    
    % Get motion information
    dx = 0;
    dy = 0;
    vel = [0 0]';
    for i = 2:size(To,3)-1
        im = TSmooth(:,:,i-1:i+1); %use the smoothed version
        
        
        [dFdx dFdy] = gradient(im(:,:,2));
        dFdt = (im(:,:,3)-im(:,:,1))/2;
        H = [dFdx(:) dFdy(:)];
        vel(:,i) = inv(H'*H)*H'*dFdt(:); %units: pixels/frame
        
        if vel(1,i) > 2 || vel(2,i) > 2
            vel(1,i) = 0; vel(2,i) = 0;
        end
        
        dx(i) = dx(i-1) + vel(1,i);
        dy(i) = dy(i-1) + vel(2,i);
        
    end
    
    dx(end+1) = dx(end);
    dy(end+1) = dy(end);
    %figure, plot(dx), hold on, plot(dy)
    
    % Shift correction
    for i = 1:size(To,3)
        To(:,:,i) = circshift_continous2(squeeze(To(:,:,i)),dx(i),dy(i));
    end
    
end
