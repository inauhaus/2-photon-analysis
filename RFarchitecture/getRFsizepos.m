function RF = getRFsizepos

global maskS cellS

%%
figure

xsize_deg = 1./getparam('s_freq2'); %width of the stimulus in degrees
ysize_deg = 1./getparam('s_freq2'); %height of stimulus in degrees
phasedom = getdomain('s_phase2');
id180 = find(phasedom == 180);

xdom = xsize_deg*phasedom/360; %convert domain from phase into degrees of visual field.
ydom = ysize_deg*phasedom/360;

Ncell = length(cellS.mukern);
for i = 1:Ncell
    
    RFx = cellS.mukern{i}(:,1);  %receptive field for ori = 0;
    RFy = cellS.mukern{i}(:,2);  %receptive field for ori = 90;
    
    %Can do the circular shift if it is 
    %RFx = circshift(RFx,[length(RFx)-id180 0]); %0 is at the 180 location, plus one more to wrap around
    RFx = flipud(RFx);  %flip it to makes left to right
    
    %RFy = circshift(RFy,[length(RFy)-id180 0]); %make 0 the middle of the screen, instead of the edge
    
    [param_x ffitx varaccx sigma] = Gaussfit(xdom,RFx',1);
    [param_y ffity varaccy sigma] = Gaussfit(ydom,RFy',1);
    
    RF.xpos(i) = param_x(1);  %xposition of RF (degrees)
    RF.ypos(i) = param_y(1); %yposition of RF (degrees)
    
    RF.xsig(i) = param_x(2); %xWidth of RF (1sigma degrees)
    RF.ysig(i) = param_y(2);  %yWidth of RF (1sigma degrees)
    
    RF.xffit(i,:) = ffitx; %store the fit
    RF.yffit(i,:) = ffity;
    
    RF.x_varacc(i) = varaccx;  %store the quality of the fit (variance accounted for)
    RF.y_varacc(i) = varaccy;
    
    figure(98)
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),i) 
    plot([RFx; RFx])
    hold on
    plot([RF.xffit(i,:) RF.xffit(i,:)])
    axis tight
    axis off
    title(num2str(round(RF.xsig(i))))
    
    figure(99)
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),i) 
    plot([RFy; RFy])
    hold on
    plot([RF.yffit(i,:) RF.yffit(i,:)])
    axis off
    title(num2str(round(RF.ysig(i))))
    axis tight

    
end

