function [Jacvec ocdomvec] = RetvsOcdom(xpos,ypos,ocdommap,D)

ocdommap = ocdommap-nanmedian(ocdommap(:));
rnge = prctile(ocdommap(:),5) - prctile(ocdommap(:),95);
ocdommap = 2*ocdommap/rnge; %normalize these units

xpos = xpos(round(D/2):D:end,round(D/2):D:end);
ypos = ypos(round(D/2):D:end,round(D/2):D:end);
ocdommap = ocdommap(round(D/2):D:end,round(D/2):D:end);

[dxdu dxdv] = gradient(xpos);
[dydu dydv] = gradient(ypos);
dydu = -dydu; dydv = -dydv;  %positive progresion is up
dxdu = dxdu(2:end-1,2:end-1); dxdv =  dxdv(2:end-1,2:end-1); %the edges aren't really derivatives
dydu = dydu(2:end-1,2:end-1); dydv = dydv(2:end-1,2:end-1);

Jac = abs(dxdu.*dydv - dxdv.*dydu);
vecX = dxdu + 1i*dxdv; vecY = dydu + 1i*dydv;
Res = abs(vecX).*exp(1i*angle(vecX)*2) + abs(vecY).*exp(1i*angle(vecY)*2);
Res = Res./(abs(vecX) + abs(vecY));
prefAxisMF = angle(Res)/2*180/pi;
id = find(prefAxisMF<0); 
prefAxisMF(id) = prefAxisMF(id)+180;
%figure,hist(prefAxisMF(:))

AR = (1-abs(Res).^2)./(1+abs(Res).^2);

dposdu = sqrt(dxdu.^2 + dydu.^2);
dposdv = sqrt(dxdv.^2 + dydv.^2);

[dODdu dODdv] = gradient(ocdommap);
dODdu = dODdu(2:end-1,2:end-1);
dODdv =  dODdv(2:end-1,2:end-1); %the edges aren't really derivatives
prefAxisOD = atan2(dODdv,dODdu)*180/pi;

%figure, hist(prefAxisOD(:))

RetODInterct = abs(angle(exp(1i*(prefAxisOD - prefAxisMF)*pi/180))*180/pi);
figure,hist(abs(RetODInterct(:))), xlabel('OcDom vs. MagFactor gradient intersection')


posvec = [dposdu dposdv];
posvec = posvec(:);
ODvec = [dODdu dODdv];
ODvec = abs(ODvec(:));

id = find(isnan(posvec.*ODvec));
posvec(id) = []; ODvec(id) = [];

figure,
scatter(ODvec,posvec,'.k')
[r p] = corrcoef(posvec,ODvec);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('|change in position| (deg )'), xlabel('|change in ocular dominance|')


%%

ocdommap = ocdommap(2:end-1,2:end-1);

ocdomvec = ocdommap(:);
Jacvec = Jac(:);

id = find(isnan(ocdomvec.*Jacvec));
ocdomvec(id) = []; Jacvec(id) = [];

figure,
scatter(Jacvec(:),ocdomvec(:),'.k')
[r p] = corrcoef(Jacvec(:),ocdomvec(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('monocularity'), xlabel('Magnification factor')

