function err = Decayfitter_handle(param)

global yy xx;

A = param(1);

B = param(2);

C = param(3);

ffit = 1./(xx*C+A) + B;

%id = find(~isnan(yy.*ffit));
%err = trimmean((ffit(id)-yy(id)).^2,20);

erry = nanmean(((ffit-yy).^2));

ffitx = (1./(yy-B) - A)/C;
errx = nanmean(((ffitx-xx).^2));

%err = min([(errx) (erry)]);
%err = errx + erry;

err = erry;
