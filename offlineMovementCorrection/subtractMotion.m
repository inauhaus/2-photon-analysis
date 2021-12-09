function To = subtractMotion(To,X)

reshapeFlag = 0;
if length(size(To)) == 3
    dim = size(To);
    To = reshape(To,[dim(1)*dim(2) dim(3)])';
    reshapeFlag = 1;
end

H = [X(:) ones(length(X),1)];
%H = ones(length(X),1)
xhat = inv(H'*H)*H'*To;

noise = H*xhat;
To = To - noise;

if reshapeFlag
    To = reshape(To',[dim(1) dim(2) size(To,1)]); %Make it 3D again if thats how it started
end