function loadMaskInfo(pathname,filename)

%This was created so that I could easily do it from outside the GUI

global G_handles maskS

S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S

maskS = S.maskS;    %f0m is a cell array with images from each condition

% set(G_handles.maskSize,'string',num2str(maskS.maskSize));
% set(G_handles.maskThresh,'string',num2str(maskS.maskThresh));
% set(G_handles.maskMorph,'string',num2str(maskS.maskMorph));
% set(G_handles.minCellArea,'string',['[' num2str(maskS.minCellArea) ']']);
