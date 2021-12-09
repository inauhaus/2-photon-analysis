function loadMaskInfo_SBX(pathname,filename)

%This was created so that I could easily do it from outside the GUI

global G_handles maskS ACQinfo

load(strcat(pathname,filename),'-mat');  %Returns the contents in the .mat under the structure S

mask = mask(:,ACQinfo.unblanked);

maskS.bwCell_sbxID{1} = mask;   %Make sure to save the original because it ID's the time courses
maskS.bwCell{1} = sign(mask);   


