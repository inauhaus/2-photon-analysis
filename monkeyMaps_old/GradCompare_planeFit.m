%%
Npc = 32;
dim = size(sfpref_raw);
Wshift = Npc/2;
[x y] = meshgrid(0:Npc-1,0:Npc-1);
H =  double([x(:) y(:) ones(Npc^2,1)]);
clear xslope_ori yslope_ori DC_ori xslope_sf yslope_sf DC_sf
for i = 1:(dim(1)/Wshift - 1)
    i
    for j = 1:(dim(2)/Wshift - 1)
        
           rows = ((i-1)*Wshift+1) : ((i-1)*Wshift+Npc);
           cols = ((j-1)*Wshift+1) : ((j-1)*Wshift+Npc);
           
           impc = sfpref_raw(rows,cols);
           y = log2(impc(:));
           
           p = inv(H'*H)*H'*y;
           
           xslope_sf(i,j) = p(1);
           yslope_sf(i,j) = p(2);
           DC_sf(i,j) = p(3);
           
           impc = oriang_raw(rows,cols);
           y = impc(:);
           
%            p = inv(H'*H)*H'*y;           
%            xslope_ori(i,j) = p(1);
%            yslope_ori(i,j) = p(2);
%            DC_ori(i,j) = p(3);                     
           
           [oriplane oripref_hat] = Planefit(H(:,1),H(:,2),double(y));
           xslope_ori(i,j) = oriplane(1);
           yslope_ori(i,j) = oriplane(2);

        
    end
end
%%
xcent = Wshift:Wshift:(dim(2)/Wshift - 1)*Wshift + .5;
ycent = Wshift:Wshift:(dim(1)/Wshift - 1)*Wshift + .5;
figure,imagesc(log2(sfpref_raw),log2([.5 8]))
hold on
quiver(xcent,ycent,xslope_ori,yslope_ori,'k')


figure,imagesc(oriang_raw),
hold on
quiver(xcent,ycent,xslope_sf,yslope_sf,'k')

sfdir = atan2(yslope_sf,xslope_sf);
oridir = atan2(yslope_ori,xslope_ori);

dphase = angle(exp(1i*sfdir)*exp(-1i*oridir))*180/pi;


bdom = -180:20:180;
bcount = histc(dphase(:),bdom); bcount = bcount(1:end-1);
b = bcount;
%b = bcount/length(idROI);

bdom_center = bdom(1:end-1)+.5*(bdom(2)-bdom(1));
figure, bar(bdom_center,b)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 1 1])
xlabel('phase of intersection')
ylabel('percentage of pixels')

%Xticks = bdom(1:2:end);
Xticks = [-180 -90 0 90 180];
set(gca,'Xtick',Xticks,'XtickLabel',Xticks,'TickDir','out'); 

b = b/sum(b);
resultant = sum(b'.*exp(1i*2*bdom_center*pi/180));
Res_angle = angle(resultant)/2*180/pi;
if Res_angle < 0
    Res_angle = 180+Res_angle;
end
Cv = 1-abs(resultant);

[pval z] = circ_rtest(bdom_center*pi/180, bcount)  %Rayleigh test of uniformity from Philip Berens' 'CircStat'
title(['Preferred intersection = ' num2str(round(Res_angle)) ';  Circ Var = ' num2str((Cv))...
    '  Rayleigh p-value = ' num2str(pval)])

figure
quiver(xcent,ycent,xslope_sf,yslope_sf,'k')
hold on
quiver(xcent,ycent,xslope_ori,yslope_ori,'r')