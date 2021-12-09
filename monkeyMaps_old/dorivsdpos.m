function [oridum prefAxisMF intsect] = dorivsdpos(xpos,ypos,oriang,D)

xpos = xpos(round(D/2):D:end,round(D/2):D:end);
ypos = ypos(round(D/2):D:end,round(D/2):D:end);
oriang = oriang(round(D/2):D:end,round(D/2):D:end);

[dxdu dxdv] = gradient(xpos);
[dydu dydv] = gradient(ypos);
dydu = -dydu; dydv = -dydv;  %positive progresion is up
dxdu = dxdu(2:end-1,2:end-1); dxdv =  dxdv(2:end-1,2:end-1); %the edges aren't really derivatives
dydu = dydu(2:end-1,2:end-1); dydv = dydv(2:end-1,2:end-1);

duposmag = sqrt(dxdu.^2 + dydu.^2);
dvposmag = sqrt(dxdv.^2 + dydv.^2);

[doridu doridv] = gradient(oriang);
doridu = doridu(2:end-1,2:end-1);
doridv = doridv(2:end-1,2:end-1); %the edges aren't really derivatives
doridu = abs(angle(exp(1i*doridu*pi/180*2))*180/pi/2); %or... 90-abs(abs(doridu)-90)
doridv = abs(angle(exp(1i*doridv*pi/180*2))*180/pi/2);

posmag = [duposmag(:); dvposmag(:)];
orimag = [doridu(:); doridv(:)];

id = find(isnan(posmag.*orimag));
posmag(id) = [];
orimag(id) = [];

figure,
scatter(orimag,posmag,'.k')
[r p] = corrcoef(posmag,orimag);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
xlabel('|change in orientation|')
ylabel('|change in position| (deg )')
axis square


%%
%get preferred axis of magnification factor in visual space
vecU = dxdu + 1i*dydu; vecV = dxdv + 1i*dydv;
%Res = abs(vecU).*exp(1i*(pi/2-angle(vecU))*2) + abs(vecV).*exp(1i*(pi/2-angle(vecV))*2);
Res = abs(vecU).*exp(1i*(pi/2+angle(vecU))*2) + abs(vecV).*exp(1i*(pi/2+angle(vecV))*2);
Res = Res./(abs(vecU) + abs(vecV));
prefAxisMF = angle(Res)/2*180/pi;
id = find(prefAxisMF<0); 
prefAxisMF(id) = prefAxisMF(id)+180;

%get preferred axis of magnification factor on the cortex
% vecX = dxdu + 1i*dxdv; vecY = dydu + 1i*dydv;
% Res = abs(vecX).*exp(1i*angle(vecX)*2) + abs(vecY).*exp(1i*angle(vecY)*2);
% Res = Res./(abs(vecX) + abs(vecY));
% prefAxisMF = angle(Res)/2*180/pi;
% id = find(prefAxisMF<0); 
% prefAxisMF(id) = prefAxisMF(id)+180;

oridum = oriang(2:end-1,2:end-1);

figure,
subplot(1,2,1)
scatter(oridum(:),prefAxisMF(:),'.k')
xlabel('Preferred Orientation')
ylabel('Short axis of visual space representation')
axis square
hold on
plot([0 180],[0 180])

subplot(1,2,2)
intsect = 90-abs(abs(oridum(:)-prefAxisMF(:))-90);
%intsect = abs(oridum(:)-prefAxisMF(:));
hist(intsect,6), xlim([0 90])

