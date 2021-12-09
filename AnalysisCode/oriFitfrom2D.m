function [param ffit varacc sigma] = oriFitfrom2D(f)

%input is the 2d power spectrum.

cols = length(f(1,:));
rows = length(f(:,1));
[xx yy] = meshgrid(1:cols,1:rows);
xx = xx-ceil(cols/2)-1;
yy = yy-ceil(rows/2);
yy = flipud(yy);

%Get rid of the "zero padding in Fourier domain"
% ymarg = mean(f,2);
% id = find(ymarg<.0001);
% f(id,:) = [];
% yy(id,:) = [];
% xx(id,:) = [];
% 
% xmarg = mean(f,1);
% id = find(xmarg<.0001);
% f(:,id) = [];
% xx(:,id) = [];
% yy(:,id) = [];
% 
% if size(yy,1) ~= size(yy,2)
%     
%     'ack; truncation fucked up'
%    asdf 
% end

interpFlag = 1;
if interpFlag
    
    [xxI yyI] = meshgrid(1:.2:cols,1:.2:rows);
    xxI = xxI-ceil(cols/2)-1;
    yyI = yyI-ceil(rows/2);
    yyI = flipud(yyI);
    
    f = interp2(xx,yy,f,xxI,yyI,'spline');
    xx = xxI;
    yy = yyI;
end


%Rotate the domain so that its continuous around the peak:
tt = atan2(yy,xx)*180/pi;
% tt = tt + (90-param(1));
% tt = angle(exp(1i*tt*pi/180))*180/pi;
% id = find(sign(tt) == -1);
% tt(id) = tt(id) + 180;
% tt = tt-90;

rr = sqrt(yy.^2 + xx.^2);

idCorner = find(rr>max(xx(:)));


rr(idCorner) = NaN;
tt(idCorner) = NaN;

tt(find(tt<0)) = tt(find(tt<0))+180;

clear oritc
bins = 0:18:180;
for i = 1:length(bins)-1;
    id = find(tt>bins(i) & tt <= bins(i+1));
    oritc(i) = mean(f(id));
end
dbins = bins(2)-bins(1);
oridom = bins(1:end-1) + dbins/2;

[param ffit varacc sigma] = Gaussfit(oridom,oritc,1);  param(2) = sigma;

% figure,plot(oritc)
% hold on
% plot(ffit)


