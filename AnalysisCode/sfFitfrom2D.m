function [param varacc sfpref sfsig sftc sfdom ffit domII] = sfFitfrom2D(f,pixperdeg)

%input is the 2d power spectrum.

cols = length(f(1,:));
rows = length(f(:,1));
[xx yy] = meshgrid(1:cols,1:rows);
xx = xx-ceil(cols/2)-1;
yy = yy-ceil(rows/2);
yy = flipud(yy);

pixperim = size(f,1); %do this before truncating below

%Rotate the domain so that its continuous around the peak:

%Get rid of the "zero padding in Fourier domain"
ymarg = mean(f,2);
id = find(ymarg<.0001);
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

rr = sqrt(yy.^2 + xx.^2);

%[idy idx] = find(rr == min(rr(:)));

rr = rr*pixperdeg/pixperim; %convert units to cyc/deg
%%
%bins = logspace(log10(.2),log10(max(rr(:))),15);
bins = linspace((0),(max(rr(:))),9);
%bins = logspace(log10(.2),log10(10),15);
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

%figure,plot(sfdom,sftc,'-o')

% [param ffit varacc sigma] = Gaussfit(sfdom,sftc,0);  param(2) = sigma;
% sfsig = param(2);
%  sfpref = param(1);

[param ffit domI domII pfit varacc] = analyzeLOGtuning(sftc(:)',sfdom,0.61);
sfsig = param.BW/2;
sfpref = param.pref;
%%
% [param ffit varacc ffitIsf domIsf] = DoGfit(sftc,sfdom);
%  [ma id] = max(ffitIsf);
%  sfpref = domIsf(id);
%  sfsig = NaN;

