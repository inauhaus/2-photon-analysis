function err = Gaborfitter2Drot_handle_nophase(param)

global RF errall;

imWy = length(RF(:,1));
imWx = length(RF(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(3)*pi/180) + y*sin(param(3)*pi/180);  %clockwise
yp = -x*sin(param(3)*pi/180) + y*cos(param(3)*pi/180);

ffitD1 = exp(-(yp).^2/(2*param(1).^2));
ffitD2 = exp(-(xp).^2/(2*param(2).^2));

Env = ffitD1.*ffitD2;

dim = size(xp);

ffit = Env.*cos(xp/dim(1)*2*pi*param(4));
ffit90 = Env.*sin(xp/dim(1)*2*pi*param(4));




RFw = abs(fft2(RF));
%RFw = RFw/norm(RFw(:));
%RFw = RFw-median(RFw(:));
%RFw = RFw/max(RFw(:));



%ffitw = sqrt(abs(fft2(ffit)).^2 + abs(fft2(ffit90)).^2);
ffitw = abs(fft2(ffit)) ;


%err = sum((ffitw(:)-RFw(:)).*(ffitw(:)-RFw(:)))

%err = sum((ffitw(id)-RFw(id)).*(ffitw(id)-RFw(id)));

%err = -norm(RFw(:)'*ffitw(:));

truncflag = 1;
if truncflag
    hh = fspecial('gaussian',size(RFw),3);
    RFwsmooth = ifft2(fft2(RFw).*abs(fft2(hh)));
    id = find(RFwsmooth>3);
else
    id = ones(1,length(RFw(:)));
end

err = -corrcoef(RFw(id),ffitw(id));
err = err(1,2);

%err = nanmean((ffitw(id)-RFw(id)).*(ffitw(id)-RFw(id)));
%err = corrcoef(ffit(:),RF(:));
%err = -err(1,2);

% id = find(~isnan(ffit(:).*RF(:)));
% ffit2 = ffit(id);
% RF2 = RF(id);
% % 
% err = corrcoef(ffit2(:),RF2(:));
% err = -err(1,2);

errall = [errall err];
