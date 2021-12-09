function colorVform(tcParams,pColordum,pSdum,pSlohi,pmeth)


%%
idB = find(~isnan(tcParams.sfprefLum.*tcParams.sfprefColor.*tcParams.sfCoMLum.*tcParams.sfCoMColor.*tcParams.sfhcoLum.*tcParams.sfhcoColor.*tcParams.sfBPLum.*tcParams.sfBPColor));
idx = find(pSdum(idB)>pSlohi(1) & pSdum(idB)<pSlohi(2));
%idx = find(pSAll(idB)<.2);
idB = idB(idx);

orioffset = -10;

pColordum = pColordum(idB);


idLum = find(pColordum>0 & pColordum<1/3);
idColorLum = find(pColordum>1/3 & pColordum<2/3);
idColor = find(pColordum>2/3 & pColordum<1);

%idB = find(pSAll>pSlohi(1) & pSAll<pSlohi(2));

switch pmeth %Which color axis (axes) to compute spatial parameters
    
    case 'mean'
        
        sfdumC{1} = log2(sqrt(tcParams.sfprefLum(idB).*tcParams.sfprefColor(idB)));
        sfdumC{2} = log2(sqrt(tcParams.sfCoMLum(idB).*tcParams.sfCoMColor(idB)));
        sfdumC{3} = log2(sqrt(tcParams.sfhcoLum(idB).*tcParams.sfhcoColor(idB)));
        sfdumC{4} = tcParams.sfBPLum(idB)/2 + tcParams.sfBPColor(idB)/2;
        sfdumC{5} = log2(sqrt(tcParams.sfBWLum(idB).*tcParams.sfBWColor(idB)));
        
        oridumC{1} = tcParams.omagColor(idB)/2+tcParams.omagLum(idB)/2;
        oridumC{2} = tcParams.dmagColor(idB)/2+tcParams.dmagLum(idB)/2;
        
        oprefdum = angle( exp(2*1i*(tcParams.oprefColor(idB)*pi/180)) + exp(2*1i*(tcParams.oprefLum(idB)*pi/180)) )*180/pi/2;
        oprefdum(find(oprefdum<0)) = oprefdum(find(oprefdum<0))+180;
        
    case 'Lum'
        
        sfdumC{1} = log2(tcParams.sfprefLum(idB));
        sfdumC{2} = log2(tcParams.sfCoMLum(idB));
        sfdumC{3} = log2(tcParams.sfhcoLum(idB));
        sfdumC{4} = tcParams.sfBPLum(idB);
        sfdumC{5} = log2(tcParams.sfBWLum(idB));
        
        oridumC{1} = tcParams.omagLum(idB);
        oridumC{2} = tcParams.dmagLum(idB);
        
        oprefdum = tcParams.oprefLum(idB);
        
        
    case 'Color'
        
        sfdumC{1} = log2(tcParams.sfprefColor(idB));
        sfdumC{2} = log2(tcParams.sfCoMColor(idB));
        sfdumC{3} = log2(tcParams.sfhcoColor(idB));
        sfdumC{4} = tcParams.sfBPColor(idB);
        sfdumC{5} = log2(tcParams.sfBWColor(idB));
        
        oridumC{1} = tcParams.omagColor(idB);
        oridumC{2} = tcParams.dmagColor(idB);
        
        oprefdum = tcParams.oprefColor(idB);
        
end




%%
ylab = {'SF peak','SF CoM','SF high pass cut off','SF BPfactor' 'SF bandwidth'};

figure
for i = 1:length(sfdumC) %loop each sf parameter
    
    sfdum = sfdumC{i};
    muSFs = [nanmean(sfdum(idLum)) nanmean(sfdum(idColorLum)) nanmean(sfdum(idColor))];
    SESFs = [nanstd(sfdum(idLum))/sqrt(length(idLum)) nanstd(sfdum(idColorLum))/sqrt(length(idColorLum)) nanstd(sfdum(idColor))/sqrt(length(idColor))];

    subplot(length(sfdumC),1,i)
    %scatter(log2(sfprefLum),pColordum,'.k'), ylim([0 1])
    scatter(pColordum,sfdum,'.k'), xlim([0 1])
    hold on, errorbar([1/6 1/2 5/6],muSFs,SESFs,'r')
    ylabel(ylab{i})
    axis square
    id = find(~isnan(pColordum.*sfdum) & pColordum>0 & pColordum<1);
    [r p] = corrcoef(pColordum(id),sfdum(id));
    title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
    [h1 p1] = ttest2(sfdum(idLum),sfdum(idColorLum));
    [h2 p2] = ttest2(sfdum(idColor),sfdum(idColorLum));
    [h3 p3] = ttest2(sfdum(idLum),sfdum(idColor))
    dum = [h1 h2 h3]
    for q = 1:3
        if dum(q)
            strsig(q) = '*';
        else
            strsig(q) = 'o';
        end
    end
    xlabel(['Color selectivity ' strsig])
    if i == 1 || i == 2 || i ==3 || i == 5
        ylim(log2([.0125 .2]))
        set(gca,'YTick',log2([.0125 .025 .05 .1 .2]))
    end
        
    
    set(gca,'XTick',[0 1/3 2/3 1])
    
    yL = get(gca,'Ylim');
    hold on
    plot([1/3 1/3],yL,'--','Color',[.5 .5 .5])
    plot([2/3 2/3],yL,'--','Color',[.5 .5 .5])
    


end


%%

ylab = {'Orientation selectivity' 'Direction selectivity'};

figure
for i = 1:length(oridumC)
    
    oridum = oridumC{i};
    muSFs = [nanmean(oridum(idLum)) nanmean(oridum(idColorLum)) nanmean(oridum(idColor))];
    SESFs = [nanstd(oridum(idLum))/sqrt(length(idLum)) nanstd(oridum(idColorLum))/sqrt(length(idColorLum)) nanstd(oridum(idColor))/sqrt(length(idColor))];

    subplot(length(oridumC),1,i)
    %scatter(log2(sfprefLum),pColordum,'.k'), ylim([0 1])
    scatter(pColordum,oridum,'.k'), xlim([0 1]),ylim([0 1])
    hold on, errorbar([1/6 1/2 5/6],muSFs,SESFs,'r')
   
    axis square
    id = find(~isnan(pColordum.*oridum) & pColordum>0 & pColordum<1);
    [r p] = corrcoef(pColordum(id),oridum(id));
    title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
    [h1 p1] = ttest2(oridum(idLum),oridum(idColorLum));
    [h2 p2] = ttest2(oridum(idColor),oridum(idColorLum));
    [h3 p3] = ttest2(oridum(idLum),sfdum(idColor))
    dum = [h1 h2 h3]
    for q = 1:3
        if dum(q)
            strsig(q) = '*';
        else
            strsig(q) = 'o';
        end
    end
    xlabel(['Color selectivity ' strsig])
%     if i == 1 || i == 2 || i ==3
%         ylim(log2([.0125 .2]))
%         set(gca,'YTick',log2([.0125 .025 .05 .1 .2]))
%     end
        
    
    set(gca,'XTick',[0 1/3 2/3 1])
    
    yL = get(gca,'Ylim');
    hold on
    plot([1/3 1/3],yL,'--','Color',[.5 .5 .5])
    plot([2/3 2/3],yL,'--','Color',[.5 .5 .5])
    
    ylabel(ylab{i})
        
end

%% 

oridiffHor = abs(oridiff(oprefdum*pi/180,orioffset*pi/180))*180/pi; %Difference from horizontal: projectors are rotated 90deg

%pColordum = pColorHiSFAll(idB);
%pColordum = pColorMidSFAll(idB);

figure
subplot(1,2,1)
scatter(oridiffHor,pColordum,'.k'), ylim([0 1])



id = find(~isnan(pColordum.*oridiffHor) & pColordum>0 & pColordum<1);
[r p] = corrcoef(pColordum(id),oridiffHor(id))
xlabel('Ori difference from horizontal')
ylabel('Color selectivity')


oriBins = [0:15:90];
clear muC SEC
for i = 1:length(oriBins)-1
    
    id = find(oridiffHor>oriBins(i) & oridiffHor<oriBins(i+1));
    muC(i) = nanmedian(pColordum(id));
    
    SEC(i) = nanstd(pColordum(id))/sqrt(length(id));
    
end
dori = oriBins(2)-oriBins(1);
%hold on, errorbar(oriBins(1:end-1)+dori/2, muC,SEC,'r')
hold on, plot(oriBins(1:end-1)+dori/2, muC,'-or')
title(['r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])