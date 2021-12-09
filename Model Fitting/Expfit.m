function [param ffit MSE] = Expfit(domain,f)

global RF
global pepORI


%%%search%%%
pepORI = domain;
RF = f;
param = expfitter;
%%%%%%%%%%%
param
ffit = param(1)*exp(-param(2)*(0:length(f)-1)) + param(3);

MSE = mean((ffit-f).*(ffit-f));

