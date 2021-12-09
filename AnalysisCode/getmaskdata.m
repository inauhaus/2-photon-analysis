function getmaskdata

global ACQinfo Tens Tens_var maskS f0m f0m_var Analyzer symbolInfo


%%
Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters
Nt = size(Tens{1},3); %Number of time points
Nc = length(Tens); %Number of conditions

for c = 1:length(Tens)-1   
    for s = 1:Nsym       
        dom{s} = getdomain(Analyzer.loops.conds{c}.symbol{s});
        symID(c,s) = find(dom{s} == Analyzer.loops.conds{c}.val{s});            
    end        
end

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);

%%

%Preallocate
for p = 1:length(celldom)
    kerndum{p} = zeros((Nc-1),Nt);
end

celldom(1) = [];
for p = 1:length(celldom)
    cellid{p} = find(masklabel == celldom(p));
end

for c = 1:Nc-1  %loop each condition
    Tdum = Tens{c};
    c
    for t = 1:Nt %loop each time point

        im = squeeze(Tdum(:,:,t));
        
        for p = 1:length(celldom) %loop each cell
                        
            kerndum{p}(c,t) = mean(im(cellid{p}));
                        
        end
                
    end
    
end

%% reshape

strc = '[';
for s = 1:Nsym
    strc = [strc ' ' num2str(length(dom{s}))];
end
strc = [strc ' ' num2str(size(Tens{c},3)) '];']
kerndims = str2num(strc)

for p = 1:length(celldom) %loop each cell
    
    kerndum{p} = kerndum{p}(:);  %Not sure if this is necessary for reshaping    
    kern{p} = reshape(kerndum{p},kerndims);

end
    
%%

msPerFrame = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
tdom = (0:size(kern{1},length(kerndims))-1)*msPerFrame/1000;

monitorFrate = 60;
Tp = getparam('t_period')/monitorFrate; %sec/cyc
samplesPerCycle = Tp/(msPerFrame/1000);
Ncycles = getparam('stim_time')/Tp; %cycles during stimulus time

colordom = getdomain('theta');
Sid = find(colordom == 90);
Mid = find(colordom == 0);
Isoid = find(colordom == 135);
Lumid = find(colordom == 45);

for s = 1:Nsym
    if strcmp(Analyzer.loops.conds{1}.symbol{s},'theta');
        colorid = s;
    end
end

stimtimeID = find(tdom>getparam('predelay') & tdom<(getparam('predelay')+getparam('stim_time')));

id0 = round(getparam('predelay')/(msPerFrame/1000));
Twin = (id0+1):id0+ceil(Ncycles)*round(samplesPerCycle);

for p = 1:length(celldom)
    
    oricurve = squeeze(mean(mean(mean(kern{p}(:,:,:,stimtimeID),1),3),4));
    [dum oribest] = max(oricurve);
    
        
    sfcurve = squeeze(mean(mean(mean(kern{p}(:,oribest,:,stimtimeID),1),2),4));
    [dum sfbest(p)] = max(sfcurve);
    sfwt(p) = sum(sfcurve.*[0 1 2]')/sum(sfcurve);

    
    Sresponse{p} = squeeze(kern{p}(Sid,oribest,sfbest(p),Twin));
    Mresponse{p} = squeeze(kern{p}(Mid,oribest,sfbest(p),Twin));
    
    IsoResponse{p} = squeeze(kern{p}(Isoid,oribest,sfbest(p),Twin));
    LumResponse{p} = squeeze(kern{p}(Lumid,oribest,sfbest(p),Twin));
    
    colorTC{p} = squeeze(mean(kern{p}(:,oribest,sfbest(p),Twin),4));
    
    %Sresponse{p} = mean(reshape(Sresponse{p},[floor(Ncycles) round(samplesPerCycle)]));        
    %Mresponse{p} = mean(reshape(Mresponse{p},[floor(Ncycles) round(samplesPerCycle)]));
    
    F1_S(p) = getF1(Sresponse{p},tdom(Twin)*1000);
    F1_M(p) = getF1(Mresponse{p},tdom(Twin)*1000);
    F1_Iso(p) = getF1(IsoResponse{p},tdom(Twin)*1000);
    F1_Lum(p) = getF1(LumResponse{p},tdom(Twin)*1000);
    
    
    colorness(p) = (colorTC{p}(Isoid)-colorTC{p}(Lumid))/(colorTC{p}(Isoid)+colorTC{p}(Lumid));
    SMness(p) = (colorTC{p}(Sid)-colorTC{p}(Mid))/(colorTC{p}(Sid)+colorTC{p}(Mid));
    
    
    dum = colorTC{p}(:)'*exp(1i*colordom(:)*2*pi/180) / sum(colorTC{p}(:));
    colorpref(p) = round(angle(dum)/2*180/pi);
    colorSelec(p) = abs(dum);
    
    if colorpref(p)<0
        colorpref(p) = colorpref(p)+180;
    end
    
    figure(10)
    subplot(7,6,p)
    plot(tdom(Twin),Sresponse{p},'m')
    hold on
    plot(tdom(Twin),Mresponse{p},'g')
    ylim([0 max([Sresponse{p}; Mresponse{p}])])
    
    figure(11)
    subplot(7,6,p)
    plot(tdom(Twin),IsoResponse{p},'r')
    hold on
    plot(tdom(Twin),LumResponse{p},'k')    
    ylim([0 max([IsoResponse{p}; LumResponse{p}])])


    figure(12)
    subplot(7,6,p)
    plot(colordom,colorTC{p},'.-k')

    
    title(num2str(colorpref(p)))
    
    
    
        
end



%%
            
function F1 = getF1(y,tdom)

    
global Analyzer

try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    %'Frame rate does not exist in Analyzer. Setting to 60Hz'
end
T = getparam('t_period')/STIMfrate*1000; %ms/cycle
ce = exp(1i*2*pi*tdom/T);

F1 = y(:)'*ce(:);
F1 = F1/length(ce);



%%
%     varflag = 0;
%     if ~isempty(Tens_var{1})
%         varflag = 1;
%     end
%     
%     xran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));
%     yran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));
%     nopix = length(yran)*length(xran);
%     
%     nc = getnoconditions;
%     
%     bflag = stimblank(getnoconditions); %if a blank exists in this experiment
%     Nloop = nc;
%     blank = [];
%     if bflag
%         Nloop = nc-1;
%         dum = f0m{end}(yran,xran);
%         blank = mean(dum(:));
%     end
%     
%     Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters
%     
%     %Make tuning curve
%     for i = 1:Nloop
%         
%         dum = f0m{i}(yran,xran);
%         tc(i) = mean(dum(:));
%         if varflag
%             dum = f0m_var{i+1}(yran,xran);
%             tc_var(i) = mean(dum(:));
%         end
%         
%     end
%     
%     
%     for i = 1:Nsym
%         allDom{i} = getdomain(symbolInfo.str{i});
%         dim(i) = length(allDom{i});
%     end
%     primDom = allDom{1};
%     
%     %Create Ndim kernel template
%     
%     switch Nsym
%         case 1
%             tcmat = zeros(dim(1),1);
%             tcourseArray = cell(dim(1),1);
%             tcourseArray_var = cell(dim(1),1);
%         case 2
%             tcmat = zeros(dim(1),dim(2));
%             tcourseArray = cell(dim(1),dim(2));
%             tcourseArray_var = cell(dim(1),dim(2));
%         case 3
%             tcmat = zeros(dim(1),dim(2),dim(3));
%             tcourseArray = cell(dim(1),dim(2),dim(3));
%             tcourseArray_var = cell(dim(1),dim(2),dim(3));
%     end
%     
%     
%     %Insert values at the correct location
%     for i = 1:length(tc)
%         tcoursedum = squeeze(sum(sum(Tens{i}(yran,xran,:),1),2))/nopix;
%         if varflag
%             tcoursedum_var = squeeze(sum(sum(Tens_var{i}(yran,xran,:),1),2))/nopix;
%         else
%             tcoursedum_var = zeros(size(tcoursedum));
%         end
%         vals = Analyzer.loops.conds{i}.val;
%         clear loc
%         for j = 1:Nsym
%             loc(j) = find(allDom{j} == vals{symbolInfo.ID(j)});
%         end
%         
%         switch Nsym
%             case 1
%                 tcmat(loc(1)) = tc(i);
%                 tcourseArray{loc(1)} = tcoursedum;
%                 tcourseArray_var{loc(1)} = tcoursedum_var;
%             case 2
%                 tcmat(loc(1),loc(2)) = tc(i);
%                 tcourseArray{loc(1),loc(2)} = tcoursedum;
%                 tcourseArray_var{loc(1),loc(2)} = tcoursedum_var;
%             case 3
%                 tcmat(loc(1),loc(2),loc(3)) = tc(i);
%                 tcourseArray{loc(1),loc(2),loc(3)} = tcoursedum;
%                 tcourseArray_var{loc(1),loc(2),loc(3)} = tcoursedum_var;
%         end
%     end
%     
%     legStr{1} = 'blank';
%     if Nsym == 3
%         
%         oppCollapse = symbolInfo.Collapse(2);
%         
%         switch oppCollapse
%             
%             case 1  %Take slice at maximum
%                 
%                 [v id] = find(tcmat(:) == max(tcmat(:)));
%                 zloc = ceil(id/(dim(1)*dim(2)));
%                 
%                 tcmat = squeeze(tcmat(:,:,zloc));
%                 tcourseArray = tcourseArray(:,:,zloc);   %Cell arrays can be indexed like this apparently
%                 tcourseArray_var = tcourseArray_var(:,:,zloc);   %Cell arrays can be indexed like this apparently
%                 
%             case 2  %Take mean over opposing parameters
%                 
%                 tcmat = squeeze(mean(tcmat,3));  %Take mean across last dimension
%                 
%                 for i = 1:dim(1)
%                     for j = 1:dim(2)
%                         tcourseNew{i,j} = 0;
%                         tcourseNew_var{i,j} = 0;
%                         for k = 1:dim(3)
%                             tcourseNew{i,j} = tcourseNew{i,j} + tcourseArray{i,j,k}/dim(3);
%                             tcourseNew_var{i,j} = tcourseNew_var{i,j} + tcourseArray_var{i,j,k}/dim(3);
%                         end
%                     end
%                 end
%                 tcourseArray = tcourseNew;
%                 tcourseArray_var = tcourseNew_var;
%                 clear tcourseNew tcourseNew_var
%                 
%         end
%         
%     end
%     
%     
% end