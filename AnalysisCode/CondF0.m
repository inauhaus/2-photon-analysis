function mu = CondF0(Tens,Flim)

global ACQinfo

Flim = Flim+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

mu = cell(1,length(Tens));

tdomdum = (0:(length(Tens{1}(:))-1)) * (ACQinfo.msPerLine/ACQinfo.pixelsPerLine);
tdomdum = reshape(tdomdum,size(Tens{1},2),size(Tens{1},1),size(Tens{1},3));
tdom = zeros(size(tdomdum,2),size(tdomdum,1),size(tdomdum,3));
for i = 1:size(tdom,3)
    tdom(:,:,i) = tdomdum(:,:,i)';
end
id = find(tdom > Flim(1) & tdom < Flim(2));
Flim_mask = zeros(size(tdom));
Flim_mask(id) = 1;

for i = 1:length(Tens)  %loop through each condition
    

    Flim_mask_dum = Flim_mask(1,1,:);

    idnan = find(isnan(Tens{i}(1,1,:)));
    Flim_mask_dum(idnan) = 0;
    
    mu{i} = nansum(Tens{i}.*Flim_mask,3)/sum(Flim_mask_dum);
    
end
