function globalPlane(xpos,ypos)

[xdom ydom] = meshgrid(1:size(xpos,2),size(xpos,1):-1:1);
H = [xdom(:) ydom(:) ones(length(xdom(:)),1)];
id = find(~isnan(xpos(:)));
slopesX = inv(H(id,:)'*H(id,:))*H(id,:)'*xpos(id);
dxdu = slopesX(1); dxdv = slopesX(2);
xposhat = H*slopesX; xposhat = reshape(xposhat,size(xpos,1),size(xpos,2));

H = meshgrid(1:size(ypos,2),size(xpos,1):-1:1);
H = [xdom(:) ydom(:) ones(length(xdom(:)),1)];  
id = find(~isnan(ypos(:)));
slopesY = inv(H(id,:)'*H(id,:))*H(id,:)'*ypos(id);
dydu = slopesY(1); dydv = slopesY(2);
yposhat = H*slopesY; yposhat = reshape(yposhat,size(ypos,1),size(ypos,2));

vecX = dxdu + 1i*dxdv; vecY = dydu + 1i*dydv;
Res = abs(vecX).*exp(1i*angle(vecX)*2) + abs(vecY).*exp(1i*angle(vecY)*2);
Res = Res./(abs(vecX) + abs(vecY));
AR = (1-abs(Res).^2)./(1+abs(Res).^2)
prefAxis = angle(Res)/2*180/pi

xposhat = xposhat - nanmean(xpos(:));
xpos2 = xpos-nanmean(xpos(:));
yposhat = yposhat - nanmean(ypos(:));
ypos2 = ypos-nanmean(ypos(:));
mi = prctile([xpos2(:); ypos2(:)],1);
ma = prctile([xpos2(:); ypos2(:)],99);

figure,
subplot(2,2,1), imagesc(xpos2,[mi ma]), axis square, colorbar
subplot(2,2,2), imagesc(ypos2,[mi ma]), axis square, colorbar
subplot(2,2,3), imagesc(xposhat,[mi ma]), axis square, colorbar
subplot(2,2,4), imagesc(yposhat,[mi ma]), axis square, colorbar
