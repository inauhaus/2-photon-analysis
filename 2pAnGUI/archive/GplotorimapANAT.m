function GplotorimapANAT(mag,ang)

global fh

[imanat] = Load2phImage(1,1);

%mag = imanat;

dim = size(ang);
mag = mag-min(mag(:));
mag = mag/max(mag(:));


% N = 5;
% kern = zeros(size(imanat));
% kern(1:N,1:N) = hann(N)*hann(N)';
% kern = kern/sum(kern(:));
% imanat = ifft2(abs(fft2(kern)).*fft2(imanat));
% 
% N = 21;
% kern = zeros(size(imanat));
% kern(1:N,1:N) = -hann(N)*hann(N)';
% kern = -kern/sum(kern(:));
% kern(ceil(N/2),ceil(N/2)) = kern(ceil(N/2),ceil(N/2)) + 1;
% imanat = ifft2(abs(fft2(kern)).*fft2(imanat));

imanat = imanat-min(imanat(:));
imanat = imanat/max(imanat(:));
imagesc(imanat),colormap gray


imfunc = ang;
imfunc = imfunc-min(imfunc(:));
imfunc = imfunc/max(imfunc(:));
imfunc = round(imfunc*63+1);
imanat = round(imanat*63+1);

grayid = gray;
hsvid = hsv;
for i = 1:dim(1)
    for j = 1:dim(2)
        imout(i,j,:) = mag(i,j).^2*hsvid(imfunc(i,j),:) + grayid(imanat(i,j),:);
    end
end

set(gcf,'Color',[1 1 1]);
imout = imout/max(imout(:));

image(imout)

axis image;

fh = gcf;
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global ACQinfo Tens1 Tens2 f0m1 f0m2 bcond

tdom = 0:length(Tens1{1}(1,1,:))-1;
tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;
if isfield(ACQinfo,'stimPredelay')
    predelay = ACQinfo.stimPredelay;
    trialtime = ACQinfo.stimTrialtime;    
    tdom = tdom-predelay;
end

rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;
%  
pos = get(event_obj,'Position'); %pos(1) is column dimension
W = 5;
xran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));
yran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));

tau = pos(2)*ACQinfo.msPerLine/1000;
tdom = tdom + tau;

figure(99)
bflag = 0;
k = 1;
for(i=0:length(f0m1)-1)
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks
        v = pepgetvalues;
        oridom(k) = v(1);
        dum1 = f0m1{i+1}(yran,xran);
        dum2 = f0m2{i+1}(yran,xran);
        tc1(i+1) = mean(dum1(:));
        tc2(i+1) = mean(dum2(:));
        k = k+1;
    else
        dum1 = f0m1{i+1}(yran,xran);
        dum2 = f0m2{i+1}(yran,xran);
        blank1 = mean(dum1(:));
        blank2 = mean(dum2(:));
        tc1(i+1) = NaN;
        tc2(i+1) = NaN; %need to do this in order to index the best/worst ori.    
        bflag = 1;
    end
end

[ma1 idma1] = max(tc1);
[mi1 idmi1] = min(tc1);
[ma2 idma2] = max(tc2);
[mi2 idmi2] = min(tc2);
tc1(find(isnan(tc1))) = []; %Get rid of blank index
tc2(find(isnan(tc2))) = [];

[oridom id] = sort(oridom);

subplot(2,2,1)
if bflag == 1
    plot([oridom(1) oridom(end)],[blank1 blank1],'k'), hold on
end
plot(oridom,tc1(id),'b'), hold on 
plot(oridom,tc1(id),'ob'), hold off
xlabel('orientation'), title('Chan 1'), xlim([0 360])

subplot(2,2,2)
if bflag == 1
    plot([oridom(1) oridom(end)],[blank2 blank2],'k'), legend('blank'), hold on
end
plot(oridom,tc2(id),'b'), hold on
plot(oridom,tc2(id),'ob'), hold off
xlabel('orientation'), title('Chan 2'), xlim([0 360])

nopix = length(yran)*length(xran);

subplot(2,2,3)
dum = squeeze(sum(sum(Tens1{idma1}(yran,xran,:),1),2))/nopix;
plot(tdom(1:end-1),dum(1:end-1)), hold on, plot(tdom(1:end-1),dum(1:end-1),'o')
hold on
dum = squeeze(sum(sum(Tens1{idmi1}(yran,xran,:),1),2))/nopix;
plot(tdom(1:end-1),dum(1:end-1),'r'), hold on, plot(tdom(1:end-1),dum(1:end-1),'or')
if isfield(ACQinfo,'stimPredelay')
    ylimits = get(gca,'Ylim');
    hold on, plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
end
hold off
xlabel('sec')

subplot(2,2,4)
dum = squeeze(sum(sum(Tens2{idma2}(yran,xran,:),1),2))/nopix;
plot(tdom(1:end-1),dum(1:end-1)), hold on, plot(tdom(1:end-1),dum(1:end-1),'o')
hold on
dum = squeeze(sum(sum(Tens2{idmi2}(yran,xran,:),1),2))/nopix;
plot(tdom(1:end-1),dum(1:end-1),'r'), hold on, plot(tdom(1:end-1),dum(1:end-1),'or')
if isfield(ACQinfo,'stimPredelay')
    ylimits = get(gca,'Ylim');
    hold on, plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
end
hold off
xlabel('sec')

tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Ori: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180) ' deg']};
       