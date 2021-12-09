function dorivsdBW(oriang,BW,D)

BW = BW(round(D/2):D:end,round(D/2):D:end);
oriang = oriang(round(D/2):D:end,round(D/2):D:end);

[dBWdu dBWdv] = gradient(BW);
dBWdu = dBWdu(2:end-1,2:end-1);
dBWdv = dBWdv(2:end-1,2:end-1); %the edges aren't really derivatives

dBWdir = atan2(dBWdv,dBWdu)*180/pi; %dire
dBWaxis = dBWdir;
id = find(dBWaxis<0); 
dBWaxis(id) = dBWaxis(id) + 180;  %make it an axis
dBWdir(id) = dBWdir(id) + 360; %0 to 360

dBWdu = abs(dBWdu); 
dBWdv = abs(dBWdv);

[doridu doridv] = gradient(oriang);
doridu = doridu(2:end-1,2:end-1);
doridv = doridv(2:end-1,2:end-1); %the edges aren't really derivatives
doridu = abs(angle(exp(1i*doridu*pi/180*2))*180/pi/2); %or... 90-abs(abs(doridu)-90)
doridv = abs(angle(exp(1i*doridv*pi/180*2))*180/pi/2);

BWmag = [dBWdu(:); dBWdv(:)];
orimag = [doridu(:); doridv(:)];

id = find(isnan(BWmag.*orimag));
BWmag(id) = [];
orimag(id) = [];

figure,
scatter(orimag,BWmag,'.k')
[r p] = corrcoef(orimag,BWmag);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
xlabel('|change in orientation|')
ylabel('|change in ON-OFF| (deg )')
axis square


oriang2 = oriang(2:end-1,2:end-1);
figure,
subplot(1,2,1), scatter(dBWdir(:),oriang2(:),'.k')
xlabel('direction of On-Off change'), ylabel('Orientation preference')
r = cxcorr(dBWdir(:), oriang2(:)*2,360);
title(['circ corr coef = ' num2str(r)])
axis square

subplot(1,2,2), scatter(dBWaxis(:),oriang2(:),'.k')
xlabel('axis of On-Off change'), ylabel('Orientation preference')
r = cxcorr(dBWaxis(:), oriang2(:),180);
title(['circ corr coef = ' num2str(r)])
axis square


