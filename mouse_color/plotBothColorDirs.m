function plotBothColorDirs(rst090,rst45135)

global ACQinfo idExamp

sfdom = rst090.sfdom;
sfdomII = rst090.sfdomII;
oridom = rst090.oridom;
%%
fr = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
tdom =(0:(size(rst090.tcourseOpt{1},2)-1))*fr/1000; 
figure
for i = 1:length(idExamp)
     
    
    subplot(length(idExamp),3,(i-1)*3+1)
    plot([oridom 360],[rst090.oritc{1}(idExamp(i),:) rst090.oritc{1}(idExamp(i),1)],'.-g')
    hold on
    plot([oridom 360],[rst090.oritc{2}(idExamp(i),:) rst090.oritc{2}(idExamp(i),1)],'.-b')
    hold on,
    plot([oridom 360],[rst45135.oritc{1}(idExamp(i),:) rst45135.oritc{1}(idExamp(i),1)],'.-k')
    hold on
    plot([oridom 360],[rst45135.oritc{2}(idExamp(i),:) rst45135.oritc{2}(idExamp(i),1)],'.-r')
    %hold on,
    %plot([0 0],[0 0],'.k')
    %axis off
        axis tight
    yl = ylim;
    yl(1) = min([0 yl(1)]);
    ylim(yl)
    set(gca,'XTick',[0 180 360])
    
    subplot(length(idExamp),3,(i-1)*3+2)
 
    semilogx((sfdom),rst090.sftc{1}(idExamp(i),:),'og')
    hold on
    semilogx((sfdomII),rst090.sftcfit{1}(idExamp(i),:),'g')
    hold on
    semilogx((sfdom),rst090.sftc{2}(idExamp(i),:),'ob')
    hold on
    semilogx((sfdomII),rst090.sftcfit{2}(idExamp(i),:),'b')    
    set(gca,'Xtick',sfdom)    
    
    semilogx((sfdom),rst45135.sftc{1}(idExamp(i),:),'ok')
    hold on
    semilogx((sfdomII),rst45135.sftcfit{1}(idExamp(i),:),'k')
    hold on
    semilogx((sfdom),rst45135.sftc{2}(idExamp(i),:),'or')
    hold on
    semilogx((sfdomII),rst45135.sftcfit{2}(idExamp(i),:),'r')    
    set(gca,'Xtick',sfdom)
    
    
    axis tight
    yl = ylim;
    yl(1) = min([0 yl(1)]);
    ylim(yl)
%     
    subplot(length(idExamp),3,(i-1)*3+3)
    plot(tdom,rst090.tcourseOpt{1}(idExamp(i),:),'g')
    hold on
    plot(tdom,rst090.tcourseOpt{2}(idExamp(i),:),'b')
    hold on
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
    xlabel('sec'), ylabel('dF/F')
    axis tight
    %title(num2str(sfpref(i)))
    
    plot(tdom,rst45135.tcourseOpt{1}(idExamp(i),:),'k')
    hold on
    plot(tdom,rst45135.tcourseOpt{2}(idExamp(i),:),'r')
    hold on
    plot([getparam('predelay') getparam('predelay')+getparam('stim_time')],[-.1 -.1],'k')
    xlabel('sec'), ylabel('dF/F')
    axis tight
    %title(num2str(sfpref(i)))
end


