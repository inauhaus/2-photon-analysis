function txt = removeCellCallback(src,event_obj)

global maskS ACQinfo

pos = round(event_obj.IntersectionPoint); %pos(1) is column dimension

masklabel = bwlabel(maskS.bwCell{1},4);

idCell = find(masklabel == masklabel(pos(2),pos(1)));

maskS.bwCell{1}(idCell) = 0;  %Remove cell

figure(40),
imagesc(maskS.im{1},'Buttondownfcn',@removeCellCallback), colormap gray
hold on
contour(maskS.bwCell{1},[.5 .5],'r','Buttondownfcn',@removeCellCallback)
hold off
axis image

txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))]};

    
    