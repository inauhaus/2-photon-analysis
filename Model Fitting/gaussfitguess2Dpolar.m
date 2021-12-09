
function G = gaussfitguess2Dpolar

%Double check Initial guesses

global RF;

y = 1:length(RF(:,1));
y = y-mean(y(:));
x = 1:length(RF(1,:));
x = x-mean(x(:));
[x y] = meshgrid(x,y);
y = flipud(y);
r = sqrt(y.^2 + x.^2);
theta = atan2(y,x)*180/pi;
id = find(sign(theta) == -1);
theta(id) = theta(id) + 180;


id = find(RF(:) == max(RF(:)));

sfguess = r(id(1));
thetaguess = theta(id(1));



G = [thetaguess 20 sfguess sfguess*2/pi min(RF(:)) max(RF(:))-min(RF(:))]
