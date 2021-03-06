pF0

global bw f0m f0m_var funcmap ACQinfo symbolInfo Analyzer G_handles idExamp flipeyebit SelectivityThresh

%set(G_handles.datadir,'string','e:\2p_data\')
%set(G_handles.analyzedir,'string','e:\2p_data\AnalyzerFiles\')

set(G_handles.datadir,'string','f:\2p_data\')
set(G_handles.analyzedir,'string','f:\2p_data\AnalyzerFiles\')

SelectivityThresh = .4;
%%%%%%%First, left hemisphere%%%%%%%%%%%

idExamp = [];

%% The first three experiments here combine two experiments to produce the
%% SF/ori maps:  an ORI/SF and a ORI/OcDom

%% expts: u001_007 u001_017; 
eno = 1;

idExamp = [ 170 170;  195 145; 223 123]; 

anim = 'ab8'; dpthresh = 1;   flipeyebit = 0;
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = CombineOcdomSFexpt(anim, dpthresh);


load(['F:\SFODfigs\Autocorr\ab8mats\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf);  %Make sure expt is loaded to get pixpermm

set(G_handles.Lwidth,'string','.1');
hh = makeMapFilter;
plotMapExamples(hh)  %shows ori and sf curves

close all
%%

eno = 2;

anim = 'ac1';
expt = 'u000_115';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['c:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1;  flipeyebit = 1; %left hemisphere
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = SfOcdom_WideField(dthresh);

load(['f:\SFODfigs\Autocorr\ac1mats1\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf);  %Make sure expt is loaded to get pixpermm

%close all
%%

eno = 3;

anim = 'ac1';
expt = 'u001_011';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['c:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1; flipeyebit = 0; %right hemisphere

idExamp = [42 90; 160 115;];  
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = SfOcdom_WideField(dthresh);

load(['f:\SFODfigs\Autocorr\ac1mats2\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf);  %Make sure expt is loaded to get pixpermm

set(G_handles.Lwidth,'string','.1');
hh = makeMapFilter;
plotMapExamples_oriocdom(hh) %shows ori and od curves
plotMapExamples(hh)  %shows ori and sf curves

close all
%%

eno = 4;

anim = 'ac1';
expt = 'u002_019';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

try
    load(['c:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])
catch
    load(['f:\Beta\f0images\f0_' anim '_' expt(2:end)])
end

dthresh = 1; flipeyebit = 0; %right hemisphere
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = SfOcdom_WideField(dthresh);

load(['f:\SFODfigs\Autocorr\ac1mats3\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf); %Make sure expt is loaded to get pixpermm

close all
%% expts: u000_084 u000_102
eno = 5;

idExamp = [];
anim = 'ab9'; 
dpthresh = 3; %For some reason this map pulls in crap if I don't set this higher
flipeyebit = 1;
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = CombineOcdomSFexpt(anim, dpthresh);

load(['f:\SFODfigs\Autocorr\ab9mats\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf);  %Make sure expt is loaded to get pixpermm

close all
%% expts:  u000_047 u000_056

%This is the same ROI as expt 57 (below)

% idExamp = [];
% anim = 'ac0'; dpthresh = 1;    flipeyebit = 1;
% [sfpref{3} ocdom{3} dprimemask{3} magodgrad{3}] = CombineOcdomSFexpt(anim, dpthresh);


%% The rest are from single experiments that contain all combinations of ORI/SF/OD 

%%
eno = 6;

anim = 'ac0';  %This one is really noisy, but still consistent (may want to  smooth more)
expt = 'u000_057';

%load expt
set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',expt)
Gsetdirectories
setGUIlabels

load(['c:\2ph_code\Beta\f0images\f0_' anim '_' expt(2:end)])

dthresh = 1; flipeyebit = 1;  %left hemisphere;
[sfpref{eno} ocdom{eno} dprimemask{eno} magodgrad{eno} ODSFdangROI{eno} ODSFdangROISEL{eno}] = SfOcdom_WideField(dthresh);

load(['f:\SFODfigs\Autocorr\ac0mats\'],'imod','imsf','imori')
[odacorr{eno} sfacorr{eno} dom{eno} odspec{eno} sfspec{eno} fdom{eno}] = get1Dacorr(imod,imsf);  %Make sure expt is loaded to get pixpermm

close all
%% Combine experiments

%..Into single vector of values ()
sfprefAll = []; ocdomAll = []; magodgradAll = [];
for i = 1:length(sfpref)    
    
    id = find(round(dprimemask{i}));  %Need to "round" because zeros are not quite zero
    
   %Don't normalize
%    sfprefAll = [sfprefAll; log2(sfpref{i}(id)) ]; 
%    
%    magodgradAll = [magodgradAll; (magodgrad{i}(id))];     
%    
%    oddum = ocdom{i}(find(dprimemask{i}));
%    oddum = oddum/std(oddum);
%    ocdomAll = [ocdomAll; oddum];
   
   %Normalize
   sfprefAll = [sfprefAll; log2(sfpref{i}(id)) ]; 
   %sfprefAll = [sfprefAll; log2(sfpref{i}(id))-mean(log2(sfpref{i}(id))) ];
   
   %Gradient is already normalized because maps were Z-scored before it
   %was taken
   magodgradAll = [magodgradAll; magodgrad{i}(id)];     
   
   oddum = ocdom{i}(id);
   oddum = oddum/std(oddum);
   ocdomAll = [ocdomAll; oddum];
   
end

binocAll = -abs(ocdomAll);
sfAll = sfprefAll;

%Plot scatter
ma = prctile(magodgradAll,99);
[mat xdom ydom] = smoothscatter(magodgradAll,sfAll,.001,.01,[prctile(magodgradAll,.1) prctile(magodgradAll,99.9)],[log2(.5) log2(8)]);
figure
subplot(1,2,1)
imagesc(xdom,ydom,-(mat)), axis xy, colormap gray, axis square
[r p] = corrcoef(magodgradAll,sfAll);
xlabel('magnitude of OD gradient (normalized)'), ylabel('SF preference (octaves)')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])

[mat xdom ydom] = smoothscatter(binocAll,sfAll,.001,.01,[prctile(binocAll,.1) prctile(binocAll,99.9)],[log2(.5) log2(8)]);
subplot(1,2,2),
imagesc(xdom,ydom,-(mat)), axis xy, colormap gray, axis square
[r p] = corrcoef(binocAll,sfAll);
xlabel('binocularity (normalized)'), ylabel('SF preference (octaves)')
title(['r = ' num2str(r(1,2)) ';  p = ' num2str(p(1,2))])

%Bin by OD gradient mag and plot errorbars
Nbins = 5; 
magodDist = NaN*zeros([length(magodgradAll)/Nbins  Nbins]);
sfDist = NaN*zeros([length(magodgradAll)/Nbins  Nbins]);
clear sfmean sfsig magodmean magodDist sfDist xlabs 
for i = 1:Nbins

    lowedge = prctile(magodgradAll,(i-1)*100/Nbins);
    highedge = prctile(magodgradAll,i*100/Nbins);

    id = find(magodgradAll > lowedge & magodgradAll <= highedge);
    
    magodDist(1:length(id),i) = magodgradAll(id);
    sfDist(1:length(id),i) = sfAll(id);
    
    magodmean(i) = mean(magodgradAll(id));
    xlabs{i} = num2str(round(magodmean(i)*1000)/1000);
    
end

subplot(1,2,1)
hold on
boxplot(sfDist,'widths',.01,'positions',magodmean,'labels',xlabs,'notch','on','symbol','')
set(gca,'XTick',magodmean,'XTicklabel',xlabs)
xlim([0 0.25])

%plot(magodmean,sfmean,'-or') %Can't even see the error bars anyway

%Bin by binocularit OD gradient mag and plot errorbars
 
clear sfmean sfsig magbinocmean xlabs sfDist
for i = 1:Nbins    
    lowedge = prctile(binocAll,(i-1)*100/Nbins);
    highedge = prctile(binocAll,i*100/Nbins);

    id = find(binocAll > lowedge & binocAll <= highedge);
    
    binocDist(1:length(id),i) = binocAll(id);
    sfDist(1:length(id),i) = sfAll(id);
    
    magbinocmean(i) = mean(binocAll(id));
    xlabs{i} = num2str(round(magbinocmean(i)*1000)/1000);
end
subplot(1,2,2)
hold on 
boxplot(sfDist,'widths',.1,'positions',magbinocmean,'labels',xlabs,'notch','on')
set(gca,'XTick',magbinocmean,'XTicklabel',xlabs)

xlim([-2.5 0])

% hold on
% errorbar(magbinocmean,sfmean,sfsig,'r')
%plot(magbinocmean,sfmean,'-or')

%Now plot histogram of intersection angles
%% Intersection histograms
intangleROIAll = []; intangleROISELAll = []; TotalA = 0;
for i = 1:length(ODSFdangROI)    
   intangleROIAll = [intangleROIAll; ODSFdangROI{i}];    
   intangleROISELAll = [intangleROISELAll; ODSFdangROISEL{i}];  
   
   TotalA = length(find(dprimemask{i})) + TotalA;
end

%length(intangleROISELAll)/TotalA   %Percentage of ROI area that had
%"reliable gradients" in both OD and SF maps

plotGradIntersection(intangleROIAll,intangleROISELAll)

%%
CCbit = 1;
idom = [1 2 3 4];
if CCbit
    odacorrSum = 0;
    sfacorrSum = 0;
    odspecSum = 0;
    sfspecSum = 0;
    figure
   for q = 1:length(idom)
       i = idom(q);
       minrange = find(dom{i}<600 & dom{i}>0);
      
       subplot(length(sfacorr),1,i)
       plot(dom{i},sfacorr{i},'r')
       %xlim([0 850]), ylim([-.6 1])
       xlabel('um')
       hold on       
       plot(dom{i},odacorr{i})
       %xlim([-850 850]), ylim([-.8 1.1])
       
       [dum sfmin(i)] = min(sfacorr{i}(minrange));
       sfmin(i) = dom{i}(minrange(sfmin(i)));
       
       [dum odmin(i)] = min(odacorr{i}(minrange));
       odmin(i) = dom{i}(minrange(odmin(i)));
       
       title(['OD min @ ' num2str(odmin(i)) '; SF min @ ' num2str(sfmin(i))])
       
       odacorrdum = interp1(dom{i},odacorr{i},dom{1});
       sfacorrdum = interp1(dom{i},sfacorr{i},dom{1});
       
       odacorrSum = odacorrdum/length(sfacorr)+odacorrSum;
       sfacorrSum = sfacorrdum/length(sfacorr)+sfacorrSum;
       
       
       odspecdum = interp1(fdom{i},odspec{i},fdom{1});
       sfspecdum = interp1(fdom{i},sfspec{i},fdom{1});
       
       odspecSum = odspecdum/length(sfspec)+odspecSum;
       sfspecSum = sfspecdum/length(sfspec)+sfspecSum;
       
   end    
   figure
   subplot(1,2,1)
   plot(dom{1},sfacorrSum/max(sfacorrSum),'r')
   xlim([0 850]), ylim([-.6 1])
   xlabel('um')
   hold on
   plot(dom{1},odacorrSum/max(odacorrSum))
   xlim([-850 850]), ylim([-.6 1.1])
   
   
   subplot(1,2,2) 
   plot(fdom{1},sfspecSum,'r')
   %xlim([0 850]), ylim([-.6 1])
   xlabel('cycles/um')
   hold on
   plot(fdom{1},odspecSum)
   xlim([0 .003])
   %xlim([-850 850]), ylim([-.6 1.1])
   
   [dum id] = max(odspecSum);
odperiod = 1/fdom{1}(id);
[dum id] = max(sfspecSum);
sfperiod = 1/fdom{1}(id);

title(['OD period = ' num2str(round(odperiod)) '; SF period = ' num2str(round(sfperiod)) ])
xlabel('cycles per micron')


   odmin./sfmin
    
end