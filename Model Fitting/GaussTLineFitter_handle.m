function err = GaussTLineFitter_handle(param)

global RF dom

xc = param(1);
sx = param(2);

A = param(3);

DC = param(4);


xx = dom-xc;

img = exp(-xx.^2./(2*sx^2));
img = A*img.*(dom+DC);

err = sum((img-RF).^2);

% err = (img-RF).^2;
% err = sort(err);
% err = sum(err(1:round(length(err)*.9)));
