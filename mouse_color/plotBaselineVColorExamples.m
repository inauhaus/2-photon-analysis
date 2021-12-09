function plotBaselineVColorExamples(rst,IsomRate,idExamp)

global ACQinfo Analyzer

basedom = logspace(log10(1/32),log10(1),6);

fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
%fr = ACQinfo.msPerLine*800;
tdom =(0:(size(rst{1}.kernT,3)-1))*fr/1000;

%%
bLHigh = 6;
bLLow = 1;
Mvec = [0.4 .5625 0];
Svec = [0.4 0 1];
%%

figure
for i = 1:length(idExamp)
    
    for bL = 1:length(rst)        
        bLcurve(bL) = rst{bL}.pM(idExamp(i));
    end    
    
    %plot baseline curve
    subplot(length(idExamp),5,(i-1)*5 + 1)
    plot(basedom,bLcurve,'o-k')
    hold on
    plot([basedom(1) basedom(end)],[0 0],'--k')
    hold on
    plot([basedom(1) basedom(end)],[1 1],'--k')
    ylim([-.2 1.4])
    
    ylabel('%M')
    if i == length(idExamp)
        xlabel('baseline')
    else i ~= length(idExamp)
        set(gca,'XTickLabel',[])
    end    
    
   %plot time courses low baseline
    subplot(length(idExamp),5,(i-1)*5 + 2)
    dum = rst{bLLow}.kernT(:,:,:,idExamp(i));    
    dum = squeeze(mean(dum,2))';  
    plot(tdom,dum(:,1),'Color',Mvec)
    hold on
    plot(tdom,dum(:,2),'Color',Svec)
    axis tight
    xlim([tdom(1) tdom(end)])
    if i == length(idExamp)
        xlabel('sec')
        set(gca,'XTick',[getparam('predelay') getparam('predelay')+getparam('stim_time')])
        set(gca,'XTickLabel',{'0';num2str(getparam('stim_time'))})
    elseif i == 1
        title('Lowest baseline')
        set(gca,'XTickLabel',[])
        legend('M response','S response')
    else
        set(gca,'XTickLabel',[])
    end
    ylabel('dF/F')
    ylimA = get(gca,'ylim');
    
    %plot time courses high baseline
    subplot(length(idExamp),5,(i-1)*5 + 3)
    dum = rst{bLHigh}.kernT(:,:,:,idExamp(i));
    dum = squeeze(mean(dum,2))';
    plot(tdom,dum(:,1),'Color',Mvec)
    hold on
    plot(tdom,dum(:,2),'Color',Svec)
    hold on
    axis tight
    %hold on
    %plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
    xlim([tdom(1) tdom(end)])
    if i == length(idExamp)
        xlabel('sec')
        set(gca,'XTick',[getparam('predelay') getparam('predelay')+getparam('stim_time')])
        set(gca,'XTickLabel',{'0';num2str(getparam('stim_time'))})
    elseif i == 1
        title('Highest baseline')
        set(gca,'XTickLabel',[])
    else 
        set(gca,'XTickLabel',[])
    end
    ylimB = get(gca,'ylim');
    
    %Get the ylims
    ylimAB(1) = min([ylimA(1) ylimB(1)])-.02;
    ylimAB(2) = max([ylimA(2) ylimB(2)])+.02;
    
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[ylimAB(1) ylimAB(1)]+.02,'k')
        
    %Set the ylims to be the same
    ylim(ylimAB)
    subplot(length(idExamp),5,(i-1)*5 + 2)
    hold on
    ylim(ylimAB)
    
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[ylimAB(1) ylimAB(1)]+.02,'k')

    
    %plot ori curve low baseline
    subplot(length(idExamp),5,(i-1)*5 + 4)
    dum = squeeze(rst{bLLow}.oricolormat(:,:,idExamp(i)))';
    dum_sig = squeeze(rst{bLLow}.oricolormat_sig(:,:,idExamp(i)))';
    dom = [rst{1}.oridom 360];     
    errorbar(dom,[dum(:,1); dum(1,1)],[dum_sig(:,1); dum_sig(1,1)],'.-','Color',Mvec)
    hold on
    errorbar(dom,[dum(:,2); dum(1,2)],[dum_sig(:,2); dum_sig(1,2)],'.-','Color',Svec)
    
%     plot([rst{1}.oridom 360],[dum(:,1); dum(1,1)],'.-','Color',Mvec)
%     hold on
%     plot([rst{1}.oridom 360],[dum(:,2); dum(1,2)],'.-','Color',Svec)
    
    xlim([0 360]), 
    axis tight
    set(gca,'XTick',[0:90:360])
    if i == length(idExamp)
        xlabel('ori')
    elseif i == 1
        title('Lowest baseline')
        set(gca,'XTickLabel',[])
    else
        set(gca,'XTickLabel',[])
    end
    ylimA = get(gca,'ylim');
    
    %plot ori curves high baseline
    subplot(length(idExamp),5,(i-1)*5 + 5)
    dum = squeeze(rst{bLHigh}.oricolormat(:,:,idExamp(i)))';
    dum_sig = squeeze(rst{bLHigh}.oricolormat_sig(:,:,idExamp(i)))';
%     plot([rst{1}.oridom 360],[dum(:,1); dum(1,1)],'.-','Color',Mvec)
%     hold on
%     plot([rst{1}.oridom 360],[dum(:,2); dum(1,2)],'.-','Color',Svec)
    dom = [rst{1}.oridom 360];     
    errorbar(dom,[dum(:,1); dum(1,1)],[dum_sig(:,1); dum_sig(1,1)],'.-','Color',Mvec)
    hold on
    errorbar(dom,[dum(:,2); dum(1,2)],[dum_sig(:,2); dum_sig(1,2)],'.-','Color',Svec)
    xlim([0 360]), 
    axis tight
    set(gca,'XTick',[0:90:360])
    if i == length(idExamp)
        xlabel('ori')
    elseif i == 1
        title('Highest baseline')
        set(gca,'XTickLabel',[])
    else 
        set(gca,'XTickLabel',[])
    end
    ylimB = get(gca,'ylim');
    
    %Get the ylims
    ylimAB(1) = min([ylimA(1) ylimB(1)])-.02;
    ylimAB(2) = max([ylimA(2) ylimB(2)])+.01;
    
    %Set the ylims to be the same
    ylim(ylimAB)
    subplot(length(idExamp),5,(i-1)*5 + 4)
    ylim(ylimAB)
    
    
end

%%

%%
clear ylimAB

figure
for i = 1:length(idExamp)
    

    
    for bL = 1:length(basedom)
        %plot time courses
        subplot(length(basedom)+1,length(idExamp),(bL-1)*length(idExamp) + i)
        dum = rst{bL}.kernT(:,:,:,idExamp(i));
        dum = squeeze(mean(dum,2))';
        plot(tdom,dum(:,1),'Color',Mvec)
        hold on
        plot(tdom,dum(:,2),'Color',Svec)
        axis tight
        xlim([tdom(1) tdom(end)])
        
        set(gca,'XTick',[getparam('predelay') getparam('predelay')+getparam('stim_time')])
        set(gca,'XTickLabel',{'0';num2str(getparam('stim_time'))})
        
        
        if bL == length(basedom)
            %xlabel('sec')
        end

        if i == 1
            ylabel('dF/F')
        end
        
        ylimAB{i}(bL,:) = get(gca,'ylim');
        if isnan(dum(1))
            ylimAB{i}(bL,:) = ylimAB{i}(bL,:)*NaN;
        end
    end
    
    for bL = 1:length(basedom)
        
        subplot(length(basedom)+1,length(idExamp),(bL-1)*length(idExamp) + i)
        
        ylimdum(1) =  min(ylimAB{i}(:,1));
        ylimdum(2) =  max(ylimAB{i}(:,2));
        plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[ylimdum(1) ylimdum(1)]+.02,'k')
        ylim(ylimdum)
        
    end
    
    
    for bL = 1:length(rst)
        bLcurve(bL) = rst{bL}.pM(idExamp(i));
    end
    
    subplot(length(basedom)+1,length(idExamp),(bL-1+1)*length(idExamp) + i)
    
    plot(basedom,bLcurve,'o-k')
    hold on
    plot([basedom(1) basedom(end)],[0 0],'--k')
    hold on
    plot([basedom(1) basedom(end)],[1 1],'--k')
    ylim([-.2 1.4])
    
    if i == 1
        ylabel('%M')
    end


end


    
   