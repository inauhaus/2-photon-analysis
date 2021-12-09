function err = NLsatfitter_handle(param)

global yy xx;


%ffit = param(1)*xx;

ffit = xx-param(2);
ffit(find(ffit>param(1))) = param(1);

%ffit = k*phi(xx-T).^n;

%id = find(~isnan(yy.*ffit));
%err = trimmean((ffit(id)-yy(id)).^2,20);

err = nanmean((ffit-yy).^2);


