function [sfpref sfsig orisig varacc_sf varacc_ori oritc orifit sftc sffit] = sforiFitfrom2D_margin(f,pixperdeg,truncflag)



%input is the 2d power spectrum.

cols = length(f(1,:));
rows = length(f(:,1));
[xx yy] = meshgrid(1:cols,1:rows);
xx = xx-ceil(cols/2)-1;
yy = yy-ceil(rows/2);
yy = flipud(yy);

pixperim = size(f,1); %do this before truncating and/or interpolating below


interpFlag = 1;
if interpFlag
    
    [xxI yyI] = meshgrid(1:.1:cols,1:.1:rows);
    xxI = xxI-ceil(cols/2)-1;
    yyI = yyI-ceil(rows/2);
    yyI = flipud(yyI);
    
    f = interp2(xx,yy,f,xxI,yyI);
    xx = xxI;
    yy = yyI;
end



%Rotate the domain so that its continuous around the peak:

if truncflag
    %Get rid of the "zero padding in Fourier domain"
    ymarg = mean(f,2);
    id = find(ymarg<.0001*max());
    f(id,:) = [];
    yy(id,:) = [];
    xx(id,:) = [];
    
    xmarg = mean(f,1);
    id = find(xmarg<.0001);
    f(:,id) = [];
    xx(:,id) = [];
    yy(:,id) = [];
    
    if size(yy,1) ~= size(yy,2)
        'ack; truncation fucked up'
        asdf
    end
end


rr = sqrt(yy.^2 + xx.^2);

idCorner = find(rr>max(xx(:)));

rr = rr*pixperdeg/pixperim; %convert units to cyc/deg

tt = atan2(yy,xx)*180/pi;
tt(find(tt<0)) = tt(find(tt<0))+180;


edirr(idCorner) = NaN;
tt(idCorner) = NaN;
%% Get sf tuning curve marginal
%bins = logspace(log10(.2),log10(max(rr(:))),15);
bins = linspace((0),(max(rr(:))),9);
bins = 0:.25:9;
%bins = logspace(log10(.2),log10(10),15);
sfTemp = zeros(size(rr));
clear sftc sfdom
k = 1;
for i = 1:length(bins)-1
    id = find(rr>bins(i) & rr <= bins(i+1));
    mu = mean(f(id));
   
    if ~isnan(mu)
        sfdom(k) = geomean(rr(id));
        sftc(k) = mu;
        k = k+1;

    end
        
end

[param sffit domI domII pfit varacc] = analyzeLOGtuning(sftc(:)',sfdom+.01,0.61);


sfsig = param.BW/2;
sfpref = param.pref;
varacc_sf = varacc;

%% Get ori tc marginal

%Rotate the domain so that its continuous around the peak:
clear oritc 
bins = 0:18:180;
for i = 1:length(bins)-1;
    id = find(tt>=bins(i) & tt <= bins(i+1));
    oritc(i) = mean(f(id));
end

dbins = bins(2)-bins(1);
oridom = bins(1:end-1) + dbins/2;

[param orifit varacc sigma] = Gaussfit(oridom,oritc,1);  param(2) = sigma;

orisig = param(2);
varacc_ori = varacc;

