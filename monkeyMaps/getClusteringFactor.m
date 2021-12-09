function clusteringFactor = getClusteringFactor(Q)

cFdom = linspace(1,4,100);
for idclust = 1:length(cFdom)
    F1F0hat = pi*exp(-(Q.SIsmearXsig/cFdom(idclust)*2*pi*Q.meas.sfpref(:)).^2/2); %analytical solution
    ME(idclust) = mean((F1F0hat-Q.meas.F1F0).^2);
end
[dum idmi] = min(ME);
clusteringFactor = cFdom(idmi);
%figure,plot(cFdom,ME)