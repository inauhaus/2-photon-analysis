%%
clear animS exptS exptSx 
global ACQinfo

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

%idExampAll = [6 22 26; 16 20 36; 12 28 57] ;
idExampAll = [6 12 22 26 ;5 16 20 36 ;5 12 28 57] ;



%exid = 2;
%ex = exdom(exid);

exdom = 1:3;
%exdom = 3
clear TC_op PW_op TC_ro PW_ro
%for eid = 1:length(exdom)
    
exdom = 1;
trial = 5;
for eid = 1:length(exdom)
    
    idExamp = idExampAll(exdom(eid),:);
    
    ex = exdom(eid);
    
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
    set(G_RChandles.dropTrials,'string','[]')
    
    
        %%%%%%%%%%%%%%
    tcourse = cellS.cellMat{trial}(idExampAll(exdom(eid),:)+1,:,1)';   
    hh = makeTemporalfilter;    
    clear tcoursex
    for i = 1:size(tcourse,2)
        
        tcoursedum = tcourse(1:680,i);
        acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
        med = median(tcoursedum);
        tcoursedum = processTcourse(tcoursedum',hh,1,acqPeriod);
        
        tcoursedum = (tcoursedum) ./ med
        
       
        tcoursex(:,i) = tcoursedum+2*i;
                
    end
    tdom = (0:length(tcoursex)-1)*acqPeriod;
    figure,plot(tdom/1000,tcoursex,'k')
    %%%%%%%%%%%%%%
    
    
    %%%%%%load randori traces %%%%%%%%
    traceroot = 'C:\CellTraces\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    load(tracepath,'cellS')
    set(G_handles.datadir,'string',[dataRoot anim '\' anim '_' expt(2:end)])
    set(G_handles.analyzedir,'string',[anaRoot anim '\' anim '_' expt(1:8)])
    Gsetdirectories
    set(G_RChandles.dropTrials,'string','[]')
    if strcmp(expt,'u009_009') & strcmp(anim,'ax3')
        set(G_RChandles.dropTrials,'string','[21:65]');
    end

    tcourse = cellS.cellMat{trial}(idExampAll(exdom(eid),:)+1,:,1)';   
    clear tcoursex
    hh = makeTemporalfilter;    
    for i = 1:size(tcourse,2)
        
        tcoursedum = tcourse(:,i);
        med = median(tcoursedum);
        acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
        tcoursedum = processTcourse(tcoursedum',hh,1,acqPeriod);
        tcoursedum = tcoursedum/med;
        
        
        tcoursex(:,i) = tcoursedum+2*i;
                
    end
    tdom = (0:length(tcoursex)-1)*acqPeriod;
    figure,plot(tdom/1000,tcoursex,'k')
    %%%%%%%%%%%%%%
    
    
    
end

%%

exp = 1;
tcourse = cellS.cellMat{1}(idExampAll(exp,:)+1,:,1)';

hh = makeTemporalfilter;

for i = 1:size(tcourse,2)
    
   tcoursedum = tcourse(:,i);
   tcoursedum
   
               tcoursedum = processTcourse(tcoursedum,hh,1,acqPeriod);              
%             id = find(tcourse<prctile(tcourse,50) & tcourse>prctile(tcourse,0));
%             tcourse = tcourse-median(tcourse); 
%             tcourse = tcourse/std(tcourse(id));
            
            noise = fftshift(xcov(tcoursedum,'unbiased'));
            noise = noise(end) - noise(end-1);
   
    
end

