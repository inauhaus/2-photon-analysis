pF0
pRev

global G_RChandles G_handles PW NB TC

set(G_RChandles.kernelLength,'string','[-200 1000]');
set(G_RChandles.LPflag,'value',0);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);
set(G_RChandles.blankNorm,'value',0);



dataRoot = 'd:\2p_data\';
anaRoot = 'c:\AnalyzerFiles\';


%%  
getAlloriposinfo

%%  Determine preferred axes of functional maps

%sfOriaxes2(dori,dBWdiff,ax,Dist,animID_d)

%%  Compare axis of distortion for each ROI to the ROIs average ori

OriprefvsRetDistortion

%%  Put all animals at each distance into a cell
% clear doriAll dsfAll doriNormAll dsfNormAll distAll
% for j = 1:length(dori{1})  %loop through each distance
%     doriAll{j} = [];
%     dposAll{j} = [];
%     doriNormAll{j} = [];
%     dposNormAll{j} = [];
%     distAll{j} = [];
%     for i = 1:length(dori) %loop through each animal
%         doriAll{j} = [doriAll{j}(:); abs(dori{i}{j}(:))];
%         dposAll{j} = [dposAll{j}(:); abs(dpos{i}{j}(:))];
%         doriNormAll{j} = [doriNormAll{j}(:); abs(doriNorm{i}{j}(:))];
%         dposNormAll{j} = [dposNormAll{j}(:); dposNorm{i}{j}(:)];
%         distAll{j} = [distAll{j}(:); dist{i}{j}(:)];
%     end
% end
%

%%

RFsize = sqrt(xsize.^2 + ysize.^2);
xypos = xpos + 1i*ypos;
BSclustering_randpos(OAng,RFsize,xypos,BWdiff,RF,abs(dori),abs(dsize),abs(dpos),abs(dBWdiff),RFEuc,Dist,animID,animID_d);  %Boot-strap clustering

%%
pairWiseOriPosPlots2((dori),dpos,(dori),dposNorm,abs(Dist),animID_d) 

pairWiseOriSizePlots((dori),sizesum,(doriNorm),sizesum,abs(Dist),animID_d)

pairWisePosSizePlots(dpos,sizesum,dposNorm,sizesum,abs(Dist),animID_d)

pairWiseOriSizePlots((dori),dsize,(doriNorm),dsize,abs(Dist),animID_d)


pairWisePosBWPlots((dpos),abs(dBWdiff),abs(dposNorm),abs(dBWdiff),abs(Dist),animID)