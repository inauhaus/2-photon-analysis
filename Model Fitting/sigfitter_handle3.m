function err = sigfitter_handle3(param)

global yy xx


xc = param(1);
C = param(2);
B = param(3);
d = abs(xx-xc);
ffit = 5*d./(C+d) + B;



err = mean((ffit-yy).^2);




