function Gplotaxismap_cellmask(mag,ang,anatflag)

global fh G_handles maskS
%mag = log(mag)
mag = mag-min(mag(:));
mag = mag/max(mag(:));

dim = size(ang);
set(gcf,'Color',[1 1 1]);

if anatflag
    
    CH = GetTrialData([1 0],1);
%     if get(G_handles.fastMotionFlag,'Value')
%         [Px_fast Py_fast] = getTrialMotion3(CH{1});
%         CH{1} = makeGeoTrx(CH{1},Px_fast,Py_fast);
%     end
    imanat = mean(CH{1}(:,:,2:end-1),3);
    %imanat = maskS.im{1};
    
    mi = prctile(imanat(:),1);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),99);
    imanat = imanat/ma;

    imfunc = ang;
    imfunc = imfunc/180;
    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);

    hsvid = hsv;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)            
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    
    x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');

else
    imfunc = ang;
    imfunc = imfunc/180;
    imfunc = round(imfunc*63+1);
    x = image(1:length(ang(1,:)),1:length(ang(:,1)),imfunc,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    
end

axis image;

fh = gcf;

colormap hsv
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global ACQinfo Fsymbol maskS cellS Analyzer G_handles symbolInfo



varflag = get(G_handles.EbarFlag,'Value');
    


nr = getnorepeats(1);

%  
pos = round(get(event_obj,'Position')); %pos(1) is column dimension

masklabel = bwlabel(maskS.bwCell{1},4);

cellID = masklabel(pos(2),pos(1)) + 1;  %starts at zero

nc = getnoconditions;
nr = getnorepeats(1);

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    nc = nc-1;
end

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters

for i = 1:Nsym
    allDom{i} = getdomain(symbolInfo.str{i});
    dim(i) = length(allDom{i});
end

if length(dim) == 1
    dim(2) = 1;
end

tcMat = zeros(dim(2),dim(1));
tcMat_sig = zeros(dim(2),dim(1));

for i = 1:nc
    miall(i) = min(min(cellS.cellMat_norm{i}(cellID,:,:)));
    maall(i) = max(max(cellS.cellMat_norm{i}(cellID,:,:)));
end
miall = min(miall);
maall = max(maall);

predelay = getparam('predelay');
trialtime = getparam('stim_time');
tdom = 0:length(cellS.muTime{i}(1,:))-1;
tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;

tdom = tdom-predelay;

if cellID ~= 1 %if not neuropil
    
    figure(99)
    
    for i = 1:nc
        for s = 1:Nsym
            idsym(s) = find(Analyzer.loops.conds{i}.val{s} == allDom{s});
        end
        if length(idsym) == 1
            idsym(2) = 1;
        end
        idplot = dim(1)*(idsym(2)-1) + idsym(1);
        if idplot == 25
            'hi'
        end
       
        
        tcMat(idsym(2),idsym(1)) = cellS.mu{i}(cellID);
        tcMat_sig(idsym(2),idsym(1)) = cellS.sig{i}(cellID);
        
            
        %subplot(5,8,idplot)
        subplot(dim(2),dim(1),idplot)
        
        if i>1
            if length(cellS.muTime{i}(1,:))>length(cellS.muTime{i-1}(1,:))
                newN = length(cellS.muTime{i-1}(1,:));
                cellS.muTime{i} = cellS.muTime{i}(:,1:newN);
                cellS.cellMat_norm{i} = cellS.cellMat_norm{i}(:,1:newN,:);
            end
        end
        
        if varflag
            fill([tdom fliplr(tdom)],[cellS.muTime{i}(cellID,:)-cellS.sigTime{i}(cellID,:) fliplr(cellS.muTime{i}(cellID,:)+cellS.sigTime{i}(cellID,:))],'r')
        else
            plot(tdom,squeeze(cellS.cellMat_norm{i}(cellID,:,:)),'r')
        end
        mi = min(min(cellS.cellMat_norm{i}(cellID,:,:)));
        hold on
        plot(tdom,cellS.muTime{i}(cellID,:),'k','LineWidth',2)
        hold on
        plot([0 trialtime],[mi mi],'k')
        ylim([miall maall])
        xlim([tdom(1) tdom(end)])
        if i == 1
            hold on
            plot([tdom(1) tdom(1)],[0 .5])
        end
            
        axis off
        hold off
            

    end
    
%     figure(98)
%     
%     if varflag
%         domplot = allDom{1}(:)*ones(1,length(tcMat(:,1)));
%         errorbar(domplot,tcMat',tcMat_sig')
%     else
%         plot(allDom{1},tcMat')
%     end
%     TCSel = abs(tcMat*exp(1i*2*allDom{1}(:)*pi/180));
%     TCSel = TCSel./sum((tcMat'))';
%     TCSel =  round(TCSel*100)/100;
%     title(['OSI = ' num2str(TCSel')])
%     
%     hold off

    
%     %plot blank responses     
%     figure(97)
%     plot(tdom,squeeze(cellS.cellMat_norm{end}(cellID,:,:)))
%     mi = min(min(cellS.cellMat_norm{end}(cellID,:,:)));
%     hold on
%     plot(tdom,cellS.muTime{end}(cellID,:),'k','LineWidth',2)
%     hold off

    
end


% figure
% for i = 1:nc
%     for j = 1:Nsym
%         val = Analyzer.loops.conds{i}.val{j};
%         idval(j) = find(allDom{j} == val);
%         tc(idval)
%     end    
%     subplot(dim(1),dim(2),i)
%     
%     tc(idval)
%     
%     if varflag
%         errorbar(tdom,cellS.mu{i}(i,:),cellS.sigTime{i}(i,:))
%     else
%         plot(tdom,cellS.mu{i}(i,:))
%     end
% end



% 
% 
% 
% 
% 
% 
% tau = pos(2)*ACQinfo.msPerLine/1000;  
% tdom = tdom + tau;
% 
% subplot(2,1,1)
% if ~isempty(blank)
%     plot([axisdom(1) axisdom(end)],[blank blank],'k'), hold on
% else 
%     %Even if no blank was shown we put a line at zero. 
%     plot([axisdom(1) axisdom(end)],[0 0],'k'), hold on  
% end
% 
% 
% if ~varflag
%     plot(axisdom,tc,'-o'), hold off
%     legend(legStr)
% else
%     errorbar(axisdom,tc(id),sqrt(tc_var(id))/SEn,'b'), hold off
% end
% xlabel(Fsymbol)
% 
% 
% %Get 'orientation selectivity index' and put into the title
% % if ~isempty(blank)
% %     tc = tc-blank;
% % end
% 
% d = size(tc);
% if d(1) == 1 || d(2) == 1
%     tcdum = tc;
% else
%     [y x] = find(tc == max(tc(:)));
%     tcdum = tc(:,x);
% end
% 
% 
% 
% TCSel = abs(tc'*exp(1i*2*axisdom'*pi/180));
% TCSel = TCSel./sum((tc))';
% TCSel =  round(TCSel*100)/100;
% title(['OSI = ' num2str(TCSel')])
% 
% Fi = 1;
% 
% subplot(2,1,2)
% if varflag
%     dum_var = squeeze(sum(sum(Tens_var{idma}(yran,xran,:),1),2))/nopix/nr;
%     errorbar(tdom(1:end-Fi),tcourseHi(1:end-Fi),sqrt(dum_var(1:end-Fi))), hold on 
% else
%     plot(tdom(1:end-Fi),tcourseHi(1:end-Fi),'-o'), hold on
% end
% 
% if varflag
%     dum_var = squeeze(sum(sum(Tens_var{idmi}(yran,xran,:),1),2))/nopix/nr;
%     errorbar(tdom(1:end-Fi),tcourseLo(1:end-Fi),sqrt(dum_var(1:end-Fi)),'r')
% else
%     plot(tdom(1:end-Fi),tcourseLo(1:end-Fi),'-or')
% end
% 
% ylimits = get(gca,'Ylim');
% plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
% 
% hold off
% xlabel('sec')
% 
% 
tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Ori: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180) ' deg']};
       