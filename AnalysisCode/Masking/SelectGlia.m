function SelectGlia

global fh

fh = gcf;

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

global maskS fh


pos = round(get(event_obj,'Position')); %pos(1) is column dimension

W2blck = 9;
W2bh = floor(W2blck/2);

%Win = hann(W2blck)*hann(W2blck)';

[x y] = meshgrid(-W2bh:W2bh,-W2bh:W2bh);
r = sqrt(x.^2 + y.^2);
sig = W2blck/3;
Win = exp(-r.^2/(2*sig^2));

blck = maskS.im{2}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
blck = blck-min(blck(:));
blck = blck/max(blck(:));
blck = blck.*Win;
idcell = find(blck>.3);

blckBW = zeros(size(blck));
blckBW(idcell) = 1;

blck = maskS.bwCell{2}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
blck = sign(blck + blckBW);

maskS.bwCell{2}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh) = blck;

figure(fh)
imagesc(maskS.im{2}), colormap gray
hold on
contour(maskS.bwCell{2},.5,'r')
hold off


txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))]};

