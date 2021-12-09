function err = expfitter4_handle(param)

global yy xx;

B = param(1);
alp = param(2);

ffit = exp(-param(2)*(xx).^param(1));

err = mean((ffit-yy).^2);

