function [RFenv RFsize RF] = getRFenvelope(tcsf,sfdom)


sfdomLin = linspace(0,sfdom(end)-.1,14); 
tcsfLin = interp1(sfdom,tcsf,sfdomLin,'spline');

Fw = [fliplr(tcsfLin(2:end)) tcsfLin(1:end-1)];
RF = real(fftshift(ifft(phi(fftshift(Fw)))));

sfdomLin = [-fliplr(sfdomLin(2:end)) sfdomLin(1:end-1)];

dsf = sfdomLin(2)-sfdomLin(1);
%%

RFquad = imag(hilbert(RF));
RFenv = sqrt(RF.^2 + RFquad.^2);
RFdom = linspace(0,2/dsf,length(RFenv)+1); RFdom = RFdom(1:end-1);

[param ffit varacc] = Gaussfit(RFdom,RFenv,0);
RFsize = 2*param(2);

% figure, plot(RFdom,RFenv)
% hold on
% plot(RFdom,RF,'r')
% hold on,
% plot(RFdom,ffit,'k')