function BWvsOcdom(BW,ocdommap,D)

ocdommap = ocdommap-nanmedian(ocdommap(:));
rnge = prctile(ocdommap(:),5) - prctile(ocdommap(:),95);
ocdommap = 2*ocdommap/rnge; %normalize these units

%Downsample
BW = BW(round(D/2):D:end,round(D/2):D:end);
ocdommap = ocdommap(round(D/2):D:end,round(D/2):D:end);

%Gradients
[dBWdu dBWdv] = gradient(BW);
dBWdu = dBWdu(2:end-1,2:end-1); 
dBWdv = dBWdv(2:end-1,2:end-1); %the edges aren't really derivatives

[dODdu dODdv] = gradient(ocdommap);
dODdu = dODdu(2:end-1,2:end-1);
dODdv =  dODdv(2:end-1,2:end-1); %the edges aren't really derivatives
prefAxisOD = atan2(dODdv,dODdu)*180/pi;

% RetODInterct = abs(angle(exp(1i*(prefAxisOD - prefAxisMF)*pi/180))*180/pi);
% figure,hist(abs(RetODInterct(:))), xlabel('OcDom vs. MagFactor gradient intersection')


BWvec = abs([dBWdu(:); dBWdv(:)]);
ODvec = abs([dODdu(:); dODdv(:)]);

id = find(isnan(BWvec.*ODvec));
BWvec(id) = []; ODvec(id) = [];

figure,
scatter(ODvec,BWvec,'.k')
[r p] = corrcoef(BWvec,ODvec);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('|ON-OFF|'), xlabel('|change in ocular dominance|')

%%

monocularity = abs(ocdommap);
OnOffdom = abs(BW);
id = find(isnan(monocularity.*OnOffdom));
monocularity(id) = []; OnOffdom(id) = [];

figure,
scatter(OnOffdom(:),monocularity(:),'.k')
[r p] = corrcoef(OnOffdom,monocularity);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
ylabel('monocularity'), xlabel('On Off dominance')

