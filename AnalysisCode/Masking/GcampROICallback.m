function txt = GcampROICallback(src,event_obj)

global maskS G_handles ACQinfo

figure(40)

[xmicperpix ymicperpix] = getImResolution;
res = geomean([xmicperpix ymicperpix]);

pos = round(event_obj.IntersectionPoint); %pos(1) is column dimension

W2blck = str2num(get(G_handles.SelectWin,'string'));
W2bh = floor(W2blck/2);


%mask = sign(phi(impctrunc-mean(impctrunc(:))));

%Win = hann(W2blck)*hann(W2blck)';

% [x y] = meshgrid(-W2bh:W2bh,-W2bh:W2bh);
% r = sqrt(x.^2 + y.^2);
% sig = W2blck/3;
% Win = exp(-r.^2/(2*sig^2));

Win = (hann(round(W2blck))*hann(round(W2blck))');
Win = fspecial('disk',W2bh);

blck = maskS.BP(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
blck = blck.*Win;

%Circular smoothing
% ir = radon(blck,0:10:170);
% h = ones(size(ir,1),1)*fspecial('gaussian',[1 size(ir,2)],2); h = h/sum(h(:));
% blck = ifft(fft(ir,[],2).*abs(fft(h,[],2)),[],2);
% blck = iradon(blck,0:10:170);
% blck = blck(1:end-1,1:end-1);

rotdom = -30:10:30; w = hann(length(rotdom)); w = w/sum(w);
blcksmooth = zeros(size(blck));
for i = 1:length(rotdom)
    blcksmooth = blcksmooth + imrotate(blck,rotdom(i),'crop')*w(i);
end
blck = blcksmooth;
blck = blck-mean(blck(:));
blck = blck/std(blck(:));

idcell = find(blck>.0);

blckBW = zeros(size(blck));
blckBW(idcell) = 1;

SE = strel('disk',1);
blckBW = imopen(blckBW,SE);

imlabel = bwlabel(blckBW);
labelid = unique(imlabel);
N = 0;
for i = 2:length(labelid)
    id = find(imlabel == labelid(i));
    N(i) = length(id);
end
[dum id] = max(N);
blckBW = zeros(size(blckBW));
blckBW(find(imlabel == labelid(id))) = 1;

blck = maskS.bwCell{1}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh);
blck = sign(blck + blckBW);



maskS.bwCell{1}(pos(2)-W2bh:pos(2)+W2bh,pos(1)-W2bh:pos(1)+W2bh) = blck;

imagesc(maskS.im{1},'Buttondownfcn',@GcampROICallback), colormap gray
hold on
contour(maskS.bwCell{1},[.5 .5],'r')
hold off
axis image

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))]};

