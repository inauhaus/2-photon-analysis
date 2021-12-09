function err = CircGaussFit_handle2(param)

global RF dom;

xc = param(1);
kappa = param(2);

A = param(3);
base = param(4);

xx = dom;

xx = xx-xc;

img = exp(kappa*cos(xx));

img = img-min(img);
img = img/max(img);

img = A*img + base;

err = sum((img-RF).^2);