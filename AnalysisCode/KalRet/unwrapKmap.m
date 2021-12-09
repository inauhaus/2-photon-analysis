function [kv kh] = unwrapKmap(kmap_vert,kmap_hor)

%%

sPerx = getparam('x_size');
sPery = getparam('y_size');

kv = kmap_vert;
kh = kmap_hor;

kvH = medfilt2(kv,[39 39]);
khH = medfilt2(kh,[39 39]);
%%%%%%%%%%%%%%%

kvdiff = kvH-kv;
id = find(kvdiff(:)>(sPerx*.5));
%kv(id) = kv(id)+sPerx;
kv(id) = sPerx/2;

khdiff = khH-kh;
id = find(khdiff(:)>(sPery*.5));
%kh(id) = kh(id)+sPery;
kh(id) = sPery/2;

%%%%%%%%%%%%%%%

kvdiff = kv-kvH;
id = find(kvdiff(:)>(sPerx*.5));
%kv(id) = kv(id)-sPerx;
kv(id) = -sPerx/2;

khdiff = kh-khH;
id = find(khdiff(:)>(sPery*.5));
%kh(id) = kh(id)-sPery;
kh(id) = -sPery/2;
