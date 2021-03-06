function plotGradIntersection(dphaseROI,dphaseROISEL)

%%

dbdom = 20;
bdom = -180:dbdom:180; 
bcount = histc(dphaseROI,bdom); 
bROI = bcount/length(dphaseROI);

bcount = histc(dphaseROISEL,bdom); 
bROISEL = bcount/length(dphaseROISEL);

figure
stairs(bdom,bROI,'k')
hold on
stairs(bdom,bROISEL,'r')

xlabel('phase of intersection')
ylabel('percentage of pixels')
%Xticks = bdom(1:2:end);
Xticks = [-180 -90 0 90 179];
set(gca,'Xtick',Xticks,'XtickLabel',Xticks,'TickDir','out'); 

bROI = bROI(1:end-1);
bROISEL = bROISEL(1:end-1);
bdom_center = bdom(1:end-1)+dbdom/2;

bROI = bROI/sum(bROI);
resultant = sum(bROI'.*exp(1i*2*bdom_center*pi/180));
Res_angleROI = angle(resultant)/2*180/pi;
%Compute intersection angle [0 90]
Res_angleROI = abs(Res_angleROI);
Res_angleROI = 90-abs(Res_angleROI-90);
%Do this instead if you want the intersection "orientation" [0 180]
%Res_angleROI = mod(Res_angleROI,180); 
SelROI = abs(resultant);

bROISEL = bROISEL/sum(bROISEL);
resultant = sum(bROISEL'.*exp(1i*2*bdom_center*pi/180));
Res_angleROISEL = angle(resultant)/2*180/pi;
%Compute intersection angle [0 90]
Res_angleROISEL = abs(Res_angleROISEL);
Res_angleROISEL = 90-abs(Res_angleROISEL-90);
%Do this instead if you want the intersection "orientation" [0 180]
%Res_angleROI = mod(Res_angleROI,180); 
SelROISEL = abs(resultant);

%[pval z] = circ_rtest(bdom_center*pi/180, bcount)  %Rayleigh test of uniformity from Philip Berens' 'CircStat'

title(['PrefROI = ' num2str(round(Res_angleROI)) ';  SelROI = ' num2str((SelROI))...
    '  PrefROISEL = ' num2str(round(Res_angleROISEL)) ';  SelROISEL = ' num2str((SelROISEL))])

xlim([-180 179])  %Hide the edge on the right side
ylim([0 max(bROISEL)+.01])
%% Now compute the "absolute angle" of intersection [0 90]

daxROI = 90-abs(abs(dphaseROI)-90); %[0 90]
daxROISEL = 90-abs(abs(dphaseROISEL)-90); %[0 90]

% daxROI = mod(dphaseROI,180); %[0 180]
% daxROI = 90-abs(daxROI-90); %[0 90]
% daxROISEL = 90-abs(abs(dphaseROISEL)-90); %[0 90]

dbdom = 10;
bdom = 0:dbdom:90;
bROI = histc(daxROI,bdom);
bROI = bROI/length(dphaseROI);

bROISEL = histc(daxROISEL,bdom);
bROISEL = bROISEL/length(dphaseROISEL);

figure, 

stairs(bdom,bROI,'k')
hold on
stairs(bdom,bROISEL,'r')
%figure, bar(bdom,b,'k')
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',[1 1 1])

xlabel('Gradient angle of intersection')
ylabel('no of pixels')

% idorth = find(daxROI > 45);
% percorthROI = length(idorth)/length(idROI);
% idorth = find(daxROISEL > 45);
% percorthROISEL = length(idorth)/length(idHist);
% title(['medROI = ' num2str(round(median(daxROI))) '; ' num2str(round(percorthROI*100)) '% > 45; '...
%     'medROISEL = ' num2str(round(median(daxROISEL))) '; ' num2str(round(percorthROISEL*100)) '% > 45'])

%[pval z] = circ_rtest(bdom_center*pi/180, bcount)  %Rayleigh test of uniformity from Philip Berens' 'CircStat'

title(['PrefROI = ' num2str(round(Res_angleROI)) ';  SelROI = ' num2str((SelROI))...
    '  PrefROISEL = ' num2str(round(Res_angleROISEL)) ';  SelROISEL = ' num2str((SelROISEL))])


Xticks = [0 45 89];
set(gca,'Xtick',Xticks,'XtickLabel',Xticks,'TickDir','out'); 

legend('ROI','ROI & SEL')
xlim([0 89]) %Hide the edge on the right
ylim([0 max(bROISEL)+.01])