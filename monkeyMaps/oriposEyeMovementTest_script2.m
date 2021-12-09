global maskS cellS TC

clear animS exptS exptSx

animS{1} = 'ax2';
animS{2} = 'ax3';
animS{3} = 'ax3';

% ID all the random orisf experiments
%exptS{1} = 'u009_009';
exptS{1} = 'u000_108'; %Keith V1
exptS{2} = 'u009_017'; %Zeus V1
exptS{3} = 'u009_009'; %Zeus V1

% ID all the randpos experiments
%exptSx{1} = 'u009_010';
exptSx{1} = 'u000_106'; %Keith V1
exptSx{2} = 'u009_015'; %Zeus V1
exptSx{3} = 'u009_010'; %Zeus V1


exdom = [3 4 7];

exdom = 2



%%

clear tdom
for i = 1:length(trialDelim)-1
    tdom(i) = trialDelim(i) * (getparam('stim_time') + getparam('predelay') + getparam('postdelay') + 1);
end

figure,plot(tdom,muX,'-ob'), hold on, plot(tdom,muY,'-or')
xlabel('seconds')
ylabel('degrees')
set(gca,'XTick',tdom)


figure,plot(muX,muY,'-ok')
xlim([-.3 .3]), ylim([-.3 .3]), axis square



std(muY)
std(muX)

range(muY)
range(muX)


%%
for exid = 1:length(exdom)

    idExamp = [];
    ex = exdom(exid);
    
    anim = animS{ex};
    expt = exptS{ex};
    exptx = exptSx{ex};
    
    %%%%%load randori mask%%%%%%%
    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    load(maskpath,'maskS')
    
    %%%%%load oripos traces and process%%%%%%%
    traceroot = 'C:\CellTraces\randpos_alignedto_randori\';
    tracepath = [traceroot anim '_' exptx(1:8) '_cellS'];
    load(tracepath,'cellS')
    
    
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' exptx(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' exptx(1:8)])
    Gsetdirectories
    
    %%%%First include all the trials%%%%%
    set(G_RChandles.dropTrials,'string','[]')
    Grandposplots3
    
    TCAll = TC;
    
    %%%%%Now loop through limited trial window%%%%%%%%
    NWin = 5;
    trialDelim = floor(linspace(1,getnotrials,NWin+1))
    clear xposMTest yposMTest TCMTest
    for tD = 1:length(trialDelim)-1 %loop each wind of trials
        
        load(maskpath,'maskS')
        load(tracepath,'cellS')
        
        tmin = trialDelim(tD);
        tmax = trialDelim(tD+1)-1;
        tDrop = 1:getnotrials;
        tDrop(tmin:tmax)=[];
        
        
        set(G_RChandles.dropTrials,'string',['[' num2str(tDrop) ']'])
        %         if strcmp(expt,'u000_092') & strcmp(anim,'ax2')
        %             set(G_RChandles.dropTrials,'string','[24:930]');
        %         end
        
        
        hh = makeTemporalfilter;
        
        trialdom = 1:1:getnotrials;
        eval(['dT = ' get(G_RChandles.dropTrials,'string')])
        trialdom(dT) = [];
        [kernels kernblank kernIm kernSig] = Ggetrandposkernel2(cellS.cellMat,trialdom,hh);
        
        Grandposplots3
        
        TCMTest{tD} = TC;
        
    end
    
    %%
    clear dY dX dP dPnorm muX muY
    
    profAvg = 0;
    for tD = 1:length(trialDelim)-1
        profAvg = profAvg+TCMTest{tD}.profileSize{1}{1}/(length(trialDelim)-1);
    end
    
   
    %now plot experimental dynamics
    for tD = 1:length(trialDelim)-1

        %muXinit = nanmean(TCMTest{1}.xpos{1}{1});
        %muYinit = nanmean(TCMTest{1}.ypos{1}{1});

        
        ypos = TCMTest{tD}.ypos{1}{2};
        xpos = TCMTest{tD}.xpos{1}{2};

        profSize = TCMTest{tD}.profileSize;
        %dY{tD} = TCMTest{tD+1}.ypos{1}{1}-TCMTest{tD}.ypos{1}{1};
        %dX{tD} = TCMTest{tD+1}.xpos{1}{1}-TCMTest{tD}.xpos{1}{1};
        %dP{tD} = sqrt(dX{tD}.^2 + dY{tD}.^2);
      
         muX(tD) = nanmedian(xpos(2:end));
         muY(tD) = nanmedian(ypos(2:end));
        
        %muX(tD) = (xpos(1));
        %muY(tD) = (ypos(1));
        
        if tD == 1
            muX_init = muX(tD);
            muY_init = muY(tD);
        end
        
        muX(tD) = muX(tD)-muX_init;
        muY(tD) = muY(tD)-muY_init;
        xpos = xpos - muX_init;
        ypos = ypos - muY_init;
            
        
        figure(103)
        subplot(length(trialDelim)-1,1,tD)
        plot([-1 1],[0 0],'--k'), hold on, plot([0 0],[-1 1],'--k')
        hold on
        plot(xpos,ypos,'.k','MarkerSize',5)
        hold on
        plot(muX(tD),muY(tD),'.r','MarkerSize',20)
        %xlim([-nanmean(profAvg) nanmean(profAvg)])
        %ylim([-nanmean(profAvg) nanmean(profAvg)])
        xlim([-.5 .5]), ylim([-.5 .5])
        axis square
        ylabel('degrees')
        title(['trials ' num2str(trialDelim(tD)) ' to ' num2str(trialDelim(tD+1)-1)...
            ' dX = ' num2str(round(muX(tD)*1000)/1000) ' dY = ' num2str(round(muY(tD)*1000)/1000)])
        

        
        
%         dPnorm{tD} = dP{tD}./profAvg;
%         hdom = [.0:.02:.15];
%         
%         figure(101)
%         subplot(length(trialDelim)-2,1,tD)
%         hist(dP{tD},hdom)
%         xlim([-hdom(2) hdom(end)])
%         title(['median = ' num2str(nanmedian(dP{tD})) 'deg'])
%         xlabel('deg')
%         
%         figure(102)
%         subplot(length(trialDelim)-2,1,tD)
%         hist(dPnorm{tD},2*hdom)
%         xlim(2*[-hdom(2) hdom(end)])
%         title(['median = ' num2str(nanmedian(dPnorm{tD})) ' (dpos/size)'])
%         xlabel('dpos/size')
    end
%%
       sqrt((muX(3)-muX(1)).^2 + (muY(3)-muY(1).^2))
        

%
end
std(muY)
std(muX)
sqrt(std(muY)^2 + std(muX)^2)
sqrt(mean((muX'-mean(muX)).^2 + (muY'-mean(muY)).^2))

range(muY)
range(muX)

