function srate = Decon(Ca)

H = [Ca(1:end-1)'];
y = Ca(2:end)';

x = inv(H'*H)*H'*y;

srate = y-H*x;

%spiketrain(find(spiketrain<50)) = 0;
%spiketrain = spiketrain-median(spiketrain);


% figure,stem(srate,'.k');
% hold on,
% plot(trace(2:end)-base,'r')
% 
% sample = srate;
% sample(find(sample==0)) = [];
% figure,hist(sample,40)