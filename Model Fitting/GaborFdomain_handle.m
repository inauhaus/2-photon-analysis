function err = GaborFdomain_handle(param)

global RF errall;

imWy = length(RF(:,1));
imWx = length(RF(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(5)*pi/180) + y*sin(param(5)*pi/180);  %clockwise
yp = -x*sin(param(5)*pi/180) + y*cos(param(5)*pi/180);

ffitD1 = exp(-(yp-param(1)).^2/(2*param(3).^2));
ffitD2 = exp(-(xp-param(2)).^2/(2*param(4).^2));
Env = ffitD1.*ffitD2;

dim = size(xp);
carr = cos(xp/dim(1)*2*pi*param(6) + param(7)*pi/180);

ffit = Env.*carr;




err = nanmean((ffit(:)-RF(:)).*(ffit(:)-RF(:)));
%err = corrcoef(ffit(:),RF(:));
%err = -err(1,2);

% id = find(~isnan(ffit(:).*RF(:)));
% ffit2 = ffit(id);
% RF2 = RF(id);
% % 
% err = corrcoef(ffit2(:),RF2(:));
% err = -err(1,2);

errall = [errall err];
