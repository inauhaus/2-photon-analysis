function err = gaussfitter_handle2(param)

global RF dom

xc = param(1);
sx = param(2);

A = param(3);

base = 0;
if length(param) == 4
    base = param(4);
end

xx = dom-xc;

img = A*exp(-xx.^2./(2*sx^2))  +  base;

err = sum((img-RF).^2);

% err = (img-RF).^2;
% err = sort(err);
% err = sum(err(1:round(length(err)*.9)));
