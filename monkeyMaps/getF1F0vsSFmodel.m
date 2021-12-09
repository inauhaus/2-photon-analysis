function [F1F0] = getF1F0vsSFmodel(SFdom,sigwx)


x = linspace(-3,3,100);
N = exp(-x.^2/(2*sigwx^2)); N = N/sum(N);
N = ones(length(SFdom),1)*N;
F1F0 = pi*abs(sum(N.*exp(1i*SFdom*x*2*pi),2));
