function [To dx dy] = TensorSubtractFlow(To)

dx = 0;
dy = 0;
dFdt_noise = zeros(size(To));
dim = size(dFdt_noise(:,:,1));
vel = [0 0]';
for i = 2:size(To,3)-1
    im = To(:,:,i-1:i+1); %use the unsmoothed version
    
    
    [dFdx dFdy] = gradient(im(:,:,2));
    dFdt = (im(:,:,3)-im(:,:,1))/2;
    H = [dFdx(:) dFdy(:)];
    vel(:,i) = inv(H'*H)*H'*dFdt(:); %pixel shift per frame
    
    dFdt_noise_dum = H*vel(:,i); %residual motion derivative
    
    dFdt_noise(:,:,i) = reshape(dFdt_noise_dum(:),dim);
    
    dx(i) = dx(i-1) + vel(1,i);
    dy(i) = dy(i-1) + vel(2,i);
    
end


%Subtract residual flow

F_noise = cumsum(dFdt_noise,3);

To = To- F_noise;


dx(end+1) = dx(end);
dy(end+1) = dy(end);
