function err = SoGfitter_handle2(param)

global RF dom

xc = param(1);
sx = param(2);
sx2 = param(5);

A = param(3);
A2 = param(4);

% base = 0;
% if length(param) == 4
%     base = param(4);
% end


img = exp(-(dom-xc).^2./(2*sx^2)) + A2*exp(-(dom).^2./(2*sx2^2)); 

img = A*img;

err = sum((img-RF).^2);

% err = (img-RF).^2;
% err = sort(err);
% err = sum(err(1:round(length(err)*.9)));
