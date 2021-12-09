function Gplotlogmap(mag,pref,anatflag,varargin)

global fh symbolInfo


if ~isempty(varargin)
    transparentMask = varargin{1};
    transparentID = find(transparentMask == 1);
end

logdom = getdomain(symbolInfo.str{1});

mi = prctile(mag(:),2);
mag = phi(mag-mi);
ma = prctile(mag(:),98);
mag = mag/ma;
mag(find(mag>1)) = 1;


% mag = mag-min(mag(:));
% mag = mag/max(mag(:));

pref = log10(pref);
dim = size(mag);
set(gcf,'Color',[1 1 1]);

id = find(isnan(pref));
mag(id) = 0;
pref(id) = min(pref(:));

if anatflag
    
    [imanat] = getExptMean([1 0],2);
    imanat = imanat{1};
    
    mi = prctile(imanat(:),1);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),99);
    imanat = imanat/ma;

    imfunc = pref;
    imfunc = imfunc-min(imfunc(:));
    imfunc = imfunc/max(imfunc(:));
    imfunc = round(imfunc*63+1);

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    
    x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');
    
    baseid = 1;

else      

%     imfunc = pref;
%     imfunc = imfunc-min(imfunc(:));
%     imfunc = imfunc/max(imfunc(:));
%     imfunc = round(imfunc*63+1);
% 
%     jetid = jet;
%     imout = zeros(dim(1),dim(2),3);
%     for i = 1:dim(1)
%         for j = 1:dim(2)
%             imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
%         end
%     end    
% 
%     if ~isempty(varargin)
%         for i = 1:3
%             imdum = imout(:,:,i);
%             imdum(transparentID) = max(imout(:))*ones(size(transparentID));
%             imout(:,:,i) = imdum;
%         end
%     end
%     
%     
%     imout = imout/max(imout(:));
%     
%     x = image(imout,'CDataMapping','direct'); 
    
    
    %imagesc(pref,'AlphaData',mag)
 
    baseid = 1;
    imfunc = pref;
    imfunc = imfunc-log10(logdom(baseid));
    imfunc = imfunc/(log10(logdom(end))-log10(logdom(baseid)));
    imfunc = round(imfunc*63+1);  %domain max at 64
    
    image(1:length(pref(1,:)),1:length(pref(:,1)),imfunc,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    colormap jet
    
end

colorbar
axis image;

fh = gcf;

logdom = logdom(baseid:end);

for i = 1:length(logdom)
    domcell{i} = logdom(i);
end
iddom = linspace(1,64,length(logdom));
colorbar('YTick',iddom,'YTickLabel',domcell)

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global ACQinfo Tens Tens_var Flim TCWin Fsymbol oppCollapse G_handles

W = TCWin;

figure(99)

varflag = get(G_handles.EbarFlag,'Value');
    
tdom = 0:length(Tens{1}(1,1,:))-1;
tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;
if isfield(ACQinfo,'stimPredelay')
    predelay = ACQinfo.stimPredelay;
    trialtime = ACQinfo.stimTrialtime;    
    tdom = tdom-predelay;
end

nr = getnorepeats(1);

SEn = sqrt(length(Flim(1):Flim(2))*nr);  %standard error normalizer for tuning curve
%  
pos = round(get(event_obj,'Position')); %pos(1) is column dimension

%%%
[tc tcourseHi tcourseLo tcourseHi_var tcourseLo_var axisdom blank legStr] = getpixeldata(pos,W);  %This does the work
%%%


tau = pos(2)*ACQinfo.msPerLine/1000;  
tdom = tdom + tau;

subplot(2,1,1)
if ~isempty(blank)
    plot([axisdom(1) axisdom(end)],[blank blank],'k'), hold on
else 
    %Even if no blank was shown we put a line at zero. 
    plot([axisdom(1) axisdom(end)],[0 0],'k'), hold on  
end

plot(axisdom,tc,'-o'), hold off
legend(legStr)

% if ~varflag
%     plot(axisdom,tc,'-o'), hold off
%     legend(legStr)
% else
%     errorbar(axisdom,tc(id),sqrt(tc_var(id))/SEn,'b'), hold off
% end
set(gca,'Xtick',round(axisdom*1000)/1000,'Xscale','log')
id = find(Fsymbol ~= '_');
xlabel(Fsymbol(id))


%Get 'orientation selectivity index' and put into the title
if ~isempty(blank)
    tc = tc-blank;
end


if oppCollapse == 3
    [y x] = find(tc == max(tc(:)));
    tcdum = tc(:,x);
else
    tcdum = tc;
end


Fi = 1;

subplot(2,1,2)
if varflag
    %dum_var = squeeze(sum(sum(Tens_var{idma}(yran,xran,:),1),2))/nopix/nr;
    dumSE = sqrt(tcourseHi_var(1:end-Fi)/nr);
    errorbar(tdom(1:end-Fi),tcourseHi(1:end-Fi),dumSE), hold on 
else
    plot(tdom(1:end-Fi),tcourseHi(1:end-Fi),'-o'), hold on
end


if varflag
    %dum_var = squeeze(sum(sum(Tens_var{idmi}(yran,xran,:),1),2))/nopix/nr;
    
    dumSE = sqrt(tcourseLo_var(1:end-Fi)/nr);
    errorbar(tdom(1:end-Fi),tcourseLo(1:end-Fi),dumSE,'r')
else
    plot(tdom(1:end-Fi),tcourseLo(1:end-Fi),'-or')
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
       ['Value: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180)]};
       

       