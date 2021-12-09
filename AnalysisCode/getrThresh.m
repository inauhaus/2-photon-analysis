function [thresh idGoodFrames] = getrThresh(rTemp,Nstd)


%Remove the linear trend over the experiment, which always is negative
H = (1:length(rTemp))';
H = [H ones(size(H))];
lparam = inv(H'*H)*H'*rTemp(:);
rTempfit = H*lparam;
rTemp2 = rTemp(:)-rTempfit(:)+lparam(2);

%Remove outliers before identifying local median
threshA = median(rTemp2)-2*std(rTemp2);
rTemp_dum = rTemp2;
rTemp_dum(find(rTemp2<threshA)) = median(rTemp2);

rTempfit_localMed = medfilt1(rTemp_dum,900);

%subtract local median from detrended rTemp.
rTemp2 = rTemp2(:)-rTempfit_localMed(:) + median(rTemp);

%

thresh = median(rTemp2)-std(rTemp2)*Nstd;

%figure(99), plot(rTemp), hold on, plot(rTemp2)
%hold on
%plot([1 length(rTemp)],[thresh thresh])

%Use statistics of the detrended rTemp in order to ID bad frames.
idGoodFrames = ones(size(rTemp));
idGoodFrames(find(rTemp2<thresh)) = 0;

