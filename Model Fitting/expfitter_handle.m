function err = expfitter_handle(param)

global RF;

dim = length(RF);

B = param(1);
alp = param(2);

A = param(3);

xx = 0:dim-1;

ffit = B*exp(-alp*xx) + A;

err = sum((ffit-RF).^2);

