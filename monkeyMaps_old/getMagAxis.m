function prefAxisMF = getMagAxis(posplanex,posplaney)


dxdu = posplanex(2);
dxdv = -posplanex(1);
dydu = posplaney(2);
dydv = -posplaney(1);


vecU = dxdu + 1i*dydu; vecV = dxdv + 1i*dydv;
%Res = abs(vecU).*exp(1i*(pi/2-angle(vecU))*2) + abs(vecV).*exp(1i*(pi/2-angle(vecV))*2);
Res = abs(vecU).*exp(1i*(pi/2+angle(vecU))*2) + abs(vecV).*exp(1i*(pi/2+angle(vecV))*2);
Res = Res./(abs(vecU) + abs(vecV));
prefAxisMF = angle(Res)/2*180/pi;
id = find(prefAxisMF<0); 
prefAxisMF(id) = prefAxisMF(id)+180;
