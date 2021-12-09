function err = sigfitter_handle2(param)

global yy xx

xc = param(1);
sig = param(2);

A = 1;
B = 0;

d = xx-xc;

ffit = A./(1 + exp(-d*sig)) + B;

err = median((ffit-yy).^2);

%%

% muY = mean(yy);
% sigY = std(yy);
% muX = mean(xx);
% sigX = std(xx);
% 
% 
% dyMat = (yy(:)-ffit(:)')/sigY;
% 
% dxMat = (xx(:)-xx(:)')/sigX;
% 
% 
% 
% r = sqrt(dxMat.^2 + dyMat.^2);
% 
% rMin = min(r,[],2); %minimum distance of each point to the fit
% 
% err = mean(rMin);





