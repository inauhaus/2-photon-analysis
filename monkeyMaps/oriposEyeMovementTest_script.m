global maskS cellS TC

clear animS exptS
animS{1} = 'ax3'; %Zeus 
animS{2} = 'ax3'; 
animS{3} = 'ax3'; 
animS{4} = 'ax3';

animS{5} = 'ax2'; %Keith
animS{6} = 'ax2';
animS{7} = 'ax2'; %Keith
animS{8} = 'ax2';

%animS{9} = 'ax3';



exptS{1} = 'u009_015'; %V1 (Zeus)
exptS{2} = 'u009_010'; %V1 (Zeus) 19 good kernels
exptS{3} = 'u009_021';  %V2 (Zeus)
exptS{4} = 'u009_022';  %V2 (Zeus) about 6 good kernels

exptS{5} = 'u000_102'; %Keith V1 (these stink)
exptS{6} = 'u000_106'; %Keith V1 about 5 kernels.  Very few cells in mask
exptS{7} = 'u000_092'; %Keith. V1 Ok.  Off center though.  On off clustering ok
exptS{8} = 'u000_095'; %Keith only one ori.  Very clear On off clustering

%exptS{9} = 'u001_035';  %V2

clear idExampAll
idExampAll{1} = [5 10 11 25];
idExampAll{2} = [10 11 12 14];
idExampAll{3} = [1 16 26 29 31];
idExampAll{4} = [15 17 18];
idExampAll{6} = [1 3 4 5 7];
idExampAll{7} = [13 15 22 28 29 42]
idExampAll{8} = []


exdom = [3 4 7];
%exdom = [2 5 6]
exdom = 2


%%
for exid = 1:length(exdom)
    
    

    try
        idExamp = idExampAll{exdom(exid)};
    catch
        idExamp = [];
    end
    
    
    ex = exdom(exid);
    
    anim = animS{ex};
    expt = exptS{ex};
    
    %load expt
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end) ''])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    
    maskroot = 'C:\CellMasks\';
    maskpath = [maskroot anim '_' expt(1:8)];
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt(1:8) '_cellS'];
    
    load(maskpath,'maskS')
    load(tracepath,'cellS')
    
    
    %%%%First include all the trials%%%%%
    set(G_RChandles.dropTrials,'string','[]')
    %         if strcmp(expt,'u000_092') & strcmp(anim,'ax2')
    %             set(G_RChandles.dropTrials,'string','[24:930]');
    %         end
    
    Grandposplots3
    
    TCAll = TC;
    
    %%%%%Now loop through limited trial window%%%%%%%%
    NWin = 3;
    trialDelim = floor(linspace(1,getnotrials,NWin+1))
    clear xposMTest yposMTest
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
        TCMTest{tD} = TC;
        
    end
    
    %
    clear dY dX dP dPnorm muX muY
    
    profAvg = 0;
    for tD = 1:length(trialDelim)-1
        profAvg = profAvg+TCMTest{tD}.profileSize{1}{1}/(length(trialDelim)-1);
    end
    
    %now plot experimental dynamics
    for tD = 1:length(trialDelim)-1
        muXinit = nanmean(TCMTest{1}.xpos{1}{1});
        muYinit = nanmean(TCMTest{1}.ypos{1}{1});

        ypos = TCMTest{tD}.ypos{1}{1}-muYinit;
        xpos = TCMTest{tD}.xpos{1}{1}-muXinit;

        profSize = TCMTest{tD}.profileSize;
        %dY{tD} = TCMTest{tD+1}.ypos{1}{1}-TCMTest{tD}.ypos{1}{1};
        %dX{tD} = TCMTest{tD+1}.xpos{1}{1}-TCMTest{tD}.xpos{1}{1};
        %dP{tD} = sqrt(dX{tD}.^2 + dY{tD}.^2);
      
        
        figure(101)
        subplot(length(trialDelim)-1,1,tD)
        plot([-1 1],[0 0],'--k'), hold on, plot([0 0],[-1 1],'--k')
        hold on
        plot(xpos,ypos,'.k','MarkerSize',5)
        hold on
        plot(nanmean(xpos),nanmean(ypos),'.r','MarkerSize',20)
        xlim([-nanmean(profAvg) nanmean(profAvg)])
        ylim([-nanmean(profAvg) nanmean(profAvg)])
        axis square
        ylabel('degrees')
        title(['trials ' num2str(trialDelim(tD)) ' to ' num2str(trialDelim(tD+1)-1)...
            ' dX = ' num2str(round(nanmean(xpos)*1000)/1000) ' dY = ' num2str(round(nanmean(ypos)*1000)/1000)])
        
        muX(tD) = nanmean(xpos);
        muY(tD) = nanmean(ypos);
        
        
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

       sqrt((muX(3)-muX(1)).^2 + (muY(3)-muY(1).^2))
        

%
end
