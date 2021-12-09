function txt = circleROICallback(src,event_obj)

global maskS G_handles

figure(40)

pos = round(event_obj.IntersectionPoint) %pos(1) is column dimension

W2blck = str2num(get(G_handles.SelectWin,'string'));
W2bh = floor(W2blck/2);

%Win = hann(W2blck)*hann(W2blck)';

[x y] = meshgrid(-W2bh:W2bh,-W2bh:W2bh);
r = sqrt(x.^2 + y.^2);
sig = W2blck/3;
Win = exp(-r.^2/(2*sig^2));

blck = zeros(size(Win));
blck(find(Win>.5)) = 1;

% blck = maskS.im{1}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
% blck = blck-min(blck(:));
% blck = blck/max(blck(:));
% blck = blck.*Win;
% idcell = find(blck>.1);
% 
% blckBW = zeros(size(blck));
% blckBW(idcell) = 1;
% blck = maskS.bwCell{1}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
% blck = sign(blck + blckBW);

dim = size(maskS.bw);

yran = pos(2)-W2bh:pos(2)+W2bh;
id = find(yran>dim(1) | yran<1);
yran(id) = [];
blck(id,:) = [];

xran = pos(1)-W2bh:pos(1)+W2bh;
id = find(xran>dim(2) | xran<1);
xran(id) = [];
blck(:,id) = [];

maskS.bwCell{1}(yran,xran) = blck;



imagesc(maskS.im{1},'Buttondownfcn',@circleROICallback), colormap gray
hold on
contour(maskS.bwCell{1},[.5 .5],'r')
hold off
axis image

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))]};

