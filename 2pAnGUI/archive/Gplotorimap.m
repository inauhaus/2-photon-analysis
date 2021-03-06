function Gplotorimap(mag,ang)

global fh
%mag = log(mag)
mag = mag-min(mag(:));
mag = mag/max(mag(:));


set(gcf,'Color',[1 1 1]);
x = image(1:length(ang(1,:)),1:length(ang(:,1)),ang*64/(180),'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');

axis image;
colormap hsv;

fh = gcf;
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global ACQinfo Tens1 Tens2 Tens1_var Tens2_var f0m1 f0m2 f0m1_var f0m2_var bcond Flim pepANA TCWin

W = TCWin;

varflag = 0;
if ~isempty(Tens1_var)
    varflag = 1;
end
    
tdom = 0:length(Tens1{1}(1,1,:))-1;
tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;
if isfield(ACQinfo,'stimPredelay')
    predelay = ACQinfo.stimPredelay;
    trialtime = ACQinfo.stimTrialtime;    
    tdom = tdom-predelay;
end

nr = pepgetnorepeats;
rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;

SEn = sqrt(length(Flim(1):Flim(2))*nr);  %standard error normalizer for tuning curve
%  
pos = round(get(event_obj,'Position')); %pos(1) is column dimension

xran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));
yran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));

tau = pos(2)*ACQinfo.msPerLine/1000;  
tdom = tdom + tau;

figure(99)
bflag = 0;
k = 1;
for i = 0:length(f0m1)-1
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks
        v = pepgetvalues;
        
        for z = 1:length(pepANA.listOfResults{i+1}.values)  %loop through each loop parameter
            if strcmp(pepANA.listOfResults{i+1}.symbols(z),'ori')
                
                paramID = z;
                oridom(i+1) = v(paramID);
                
            else
                pdumID = z;
                pdum(i+1) = v(pdumID);
                
            end
        end

        dum1 = f0m1{i+1}(yran,xran);
        dum2 = f0m2{i+1}(yran,xran);
        tc1(i+1) = mean(dum1(:));
        tc2(i+1) = mean(dum2(:));
        if varflag
            dum1 = f0m1_var{i+1}(yran,xran);
            dum2 = f0m2_var{i+1}(yran,xran);
            tc1_var(i+1) = mean(dum1(:));
            tc2_var(i+1) = mean(dum2(:));
        end
        k = k+1;
    else
        dum1 = f0m1{i+1}(yran,xran);
        dum2 = f0m2{i+1}(yran,xran);
        blank1 = mean(dum1(:));
        blank2 = mean(dum2(:));
        tc1(i+1) = NaN;
        tc2(i+1) = NaN;     %need to do this in order to index the best/worst ori.   
        if varflag
            tc1_var(i+1) = NaN;
            tc2_var(i+1) = NaN;
        end
        pdum(i+1) = NaN;
        oridom(i+1) = NaN;
        bflag = 1;
    end
end

[ma1 idma1] = max(tc1);
if length(pepANA.listOfResults{1}.values) > 1  %if multiple params in looper
    dumslice = find(pdum == pdum(idma1));
    mi1 = min(tc1(dumslice));
    idmi1 = find(tc1 == mi1);
    idmi1 = idmi1(1);
    tc1 = tc1(dumslice);
    oridom1 = oridom(dumslice);
else
    [mi1 idmi1] = min(tc1);
    tc1(find(isnan(tc1))) = []; %Get rid of blank index
    oridom1 = oridom;
    oridom1(find(isnan(oridom1))) = [];
end

[ma2 idma2] = max(tc2);
if length(pepANA.listOfResults{1}.values) > 1
    dumslice = find(pdum == pdum(idma2));
    mi2 = min(tc2(dumslice));
    idmi2 = find(tc2 == mi2);
    idmi2 = idmi2(1);
    tc2 = tc2(dumslice);
    oridom2 = oridom(dumslice);
else
    [mi2 idmi2] = min(tc2);
    tc2(find(isnan(tc2))) = [];
    oridom2 = oridom;
    oridom2(find(isnan(oridom2))) = [];
end

[oridom1 id1] = sort(oridom1);
[oridom2 id2] = sort(oridom2);

subplot(2,2,1)
if bflag == 1
    plot([oridom1(1) oridom1(end)],[blank1 blank1],'k'), hold on
end

if ~varflag
    plot(oridom1,tc1(id1),'b'), hold on 
    plot(oridom1,tc1(id1),'ob'), hold off
else
    errorbar(oridom1,tc1(id1),sqrt(tc1_var(id1))/SEn,'b'), hold off
end
xlabel('orientation'), title('Chan 1'), xlim([0 360])

subplot(2,2,2)
if bflag == 1
    plot([oridom2(1) oridom2(end)],[blank2 blank2],'k'), legend('blank'), hold on
end

if ~varflag
    plot(oridom2,tc2(id2),'b'), hold on
    plot(oridom2,tc2(id2),'ob'), hold off
else
    errorbar(oridom2,tc2(id2),sqrt(tc2_var(id2))/SEn,'b'), hold off
end
xlabel('orientation'), title('Chan 2'), xlim([0 360])

nopix = length(yran)*length(xran);
Fi = 2;

subplot(2,2,3)
dum = squeeze(sum(sum(Tens1{idma1}(yran,xran,:),1),2))/nopix;
if varflag
    dum_var = squeeze(sum(sum(Tens1_var{idma1}(yran,xran,:),1),2))/nopix/nr;
    errorbar(tdom(1:end-Fi),dum(1:end-Fi),sqrt(dum_var(1:end-Fi))), hold on 
else
    plot(tdom(1:end-Fi),dum(1:end-Fi)), hold on, plot(tdom(1:end-Fi),dum(1:end-Fi),'o')
end
dum = squeeze(sum(sum(Tens1{idmi1}(yran,xran,:),1),2))/nopix;
if varflag
    dum_var = squeeze(sum(sum(Tens1_var{idmi1}(yran,xran,:),1),2))/nopix/nr;
    errorbar(tdom(1:end-Fi),dum(1:end-Fi),sqrt(dum_var(1:end-Fi)),'r')
else
    plot(tdom(1:end-Fi),dum(1:end-Fi),'r'), plot(tdom(1:end-Fi),dum(1:end-Fi),'or')
end
if isfield(ACQinfo,'stimPredelay')
    ylimits = get(gca,'Ylim');
    plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
end
hold off
xlabel('sec')

subplot(2,2,4)
dum = squeeze(sum(sum(Tens2{idma2}(yran,xran,:),1),2))/nopix;
if varflag
    dum_var = squeeze(sum(sum(Tens2_var{idma1}(yran,xran,:),1),2))/nopix/nr;
    errorbar(tdom(1:end-Fi),dum(1:end-Fi),sqrt(dum_var(1:end-Fi))), hold on
else
    plot(tdom(1:end-Fi),dum(1:end-Fi)), hold on, plot(tdom(1:end-Fi),dum(1:end-Fi),'o')
end
dum = squeeze(sum(sum(Tens2{idmi2}(yran,xran,:),1),2))/nopix;
if varflag
    dum_var = squeeze(sum(sum(Tens2_var{idmi1}(yran,xran,:),1),2))/nopix/nr;
    errorbar(tdom(1:end-Fi),dum(1:end-Fi),sqrt(dum_var(1:end-Fi)),'r')
else
    plot(tdom(1:end-Fi),dum(1:end-Fi),'r'), plot(tdom(1:end-Fi),dum(1:end-Fi),'or')
end
if isfield(ACQinfo,'stimPredelay')
    ylimits = get(gca,'Ylim');
    plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
end
hold off
xlabel('sec')

tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Ori: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180) ' deg']};
       