function dBWvsdpos(xpos,ypos,BW,D)

xpos = xpos(round(D/2):D:end,round(D/2):D:end);
ypos = ypos(round(D/2):D:end,round(D/2):D:end);
BW = BW(round(D/2):D:end,round(D/2):D:end);

[dxdu dxdv] = gradient(xpos);
[dydu dydv] = gradient(ypos);
dxdu = dxdu(2:end-1,2:end-1); dxdv =  dxdv(2:end-1,2:end-1); %the edges aren't really derivatives
dydu = dydu(2:end-1,2:end-1); dydv = dydv(2:end-1,2:end-1);
duposmag = sqrt(dxdu.^2 + dydu.^2);
dvposmag = sqrt(dxdv.^2 + dydv.^2);

[dBWdu dBWdv] = gradient(BW);
dBWdu = dBWdu(2:end-1,2:end-1);
dBWdv = dBWdv(2:end-1,2:end-1); %the edges aren't really derivatives
dBWdu = abs(dBWdu); 
dBWdv = abs(dBWdv);

posmag = [duposmag(:); dvposmag(:)];
BWmag = [dBWdu(:); dBWdv(:)];

id = find(isnan(posmag.*BWmag));
posmag(id) = [];
BWmag(id) = [];

figure,
scatter(BWmag,posmag,'.k')
[r p] = corrcoef(posmag,BWmag);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
xlabel('|change in ON-OFF|')
ylabel('|change in position| (deg )')
axis square


