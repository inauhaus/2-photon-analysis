function Gsetdirectories

%sets data and analyzer directories based on what is in the GUI.  Its also
%loads the acquisition parameters

global G_handles twophDATADIR ACQinfo Analyzer

%Get analyzer file
anapath = get(G_handles.analyzedir,'String'); %partial path for analyzer file
anapath = [anapath '.analyzer']
load(anapath,'Analyzer','-mat')

%% Get acquisition structure and set global

datapath = get(G_handles.datadir,'String'); %partial path for .tiff files 

twophDATADIR = datapath;  %Make path for .sbx files

load(twophDATADIR,'info')

ACQinfo.SBInfo = info;
fixSBsyncs

ACQinfo.msPerLine = 1000/ACQinfo.SBInfo.resfreq;
ACQinfo.linesPerFrame = ACQinfo.SBInfo.config.lines;
%ACQinfo.pixelsPerLine = ACQinfo.SBInfo.sz(2);

ACQinfo.unblanked = 54:729;  %I just plotted and looked at where it falls off. 
ACQinfo.pixelsPerLine = length(ACQinfo.unblanked);

[xmicperpix ymicperpix] = getImResolution;
ACQinfo.xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ACQinfo.ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
