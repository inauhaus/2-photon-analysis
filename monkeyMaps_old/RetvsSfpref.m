function [Jacvec SFvec] = RetvsSfpref(xpos,ypos,sfmap,D)

xpos = xpos(round(D/2):D:end,round(D/2):D:end);
ypos = ypos(round(D/2):D:end,round(D/2):D:end);
sfmap = sfmap(round(D/2):D:end,round(D/2):D:end);

[dxdu dxdv] = gradient(xpos);
[dydu dydv] = gradient(ypos);
dydu = -dydu; dydv = -dydv;  %positive progresion is up
dxdu = dxdu(2:end-1,2:end-1); dxdv =  dxdv(2:end-1,2:end-1); %the edges aren't really derivatives
dydu = dydu(2:end-1,2:end-1); dydv = dydv(2:end-1,2:end-1);

sfmap = sfmap(2:end-1,2:end-1);  %just to make it the same size

dposdu = sqrt(dxdu.^2 + dydu.^2);
dposdv = sqrt(dxdv.^2 + dydv.^2);

posvec = [dposdu dposdv];
posvec = posvec(:);

SFvec = log2([sfmap(:); sfmap(:)]);

id = find(isnan(posvec.*SFvec));
posvec(id) = []; SFvec(id) = [];

figure,
scatter(SFvec,posvec,'.k')
[r p] = corrcoef(posvec,SFvec);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('|change in position| (deg )'), xlabel('Spatial frequency preference (octaves)')

%%

Jac = abs(dxdu.*dydv - dxdv.*dydu);

SFvec = log2(sfmap(:));
Jacvec = Jac(:);

id = find(isnan(Jacvec.*SFvec));
Jacvec(id) = []; SFvec(id) = [];

figure,
scatter(Jacvec(:),SFvec(:),'.k')
[r p] = corrcoef(Jacvec(:),SFvec(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('SF preference (octaves)'), xlabel('Magnification factor')


