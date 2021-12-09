function err = NLsat_softfitter_handle(param)

global yy xx;


ffit = param(1)*xx;

%ffit = ffit-param(2);
ffit(find(ffit>param(2))) = param(2);


ffit = smoothFit(xx,ffit,param,.1);

%ffit = k*phi(xx-T).^n;

%id = find(~isnan(yy.*ffit));
%err = trimmean((ffit(id)-yy(id)).^2,20);

%a = (ffit(:));
%b = (yy(:));

a = log(ffit(:));
b = log(yy(:));

err = nanmean((a-b).^2);



