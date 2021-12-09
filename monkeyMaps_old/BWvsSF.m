function BWvsSF(BW,sfmap,D)

%Downsample
BW = BW(round(D/2):D:end,round(D/2):D:end);
sfmap = sfmap(round(D/2):D:end,round(D/2):D:end);

%Gradients
[dBWdu dBWdv] = gradient(BW);
dBWdu = dBWdu(2:end-1,2:end-1); 
dBWdv = dBWdv(2:end-1,2:end-1); %the edges aren't really derivatives

[dsfdu dsfdv] = gradient(log2(sfmap));
dsfdu = dsfdu(2:end-1,2:end-1);
dsfdv =  dsfdv(2:end-1,2:end-1); %the edges aren't really derivatives
prefAxissf = atan2(dsfdv,dsfdu)*180/pi;

% RetsfInterct = abs(angle(exp(1i*(prefAxissf - prefAxisMF)*pi/180))*180/pi);
% figure,hist(abs(RetsfInterct(:))), xlabel('sf vs. MagFactor gradient intersection')


BWvec = abs([dBWdu(:); dBWdv(:)]);
sfvec = abs([dsfdu(:); dsfdv(:)]);

id = find(isnan(BWvec.*sfvec));
BWvec(id) = []; sfvec(id) = [];

figure,
scatter(sfvec,BWvec,'.k')
[r p] = corrcoef(BWvec,sfvec);
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
xlabel('|change in SF preference|'), ylabel('|ON-OFF|')

%%
OnOffdom = abs(BW);
id = find(isnan(sfmap.*OnOffdom));
sfmap(id) = []; OnOffdom(id) = [];

figure,
scatter(OnOffdom(:),log2(sfmap(:)),'.k')
[r p] = corrcoef(OnOffdom,log2(sfmap));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
xlabel('On Off dominance'), ylabel('sf preference (octaves)')

