function varargout = pF0(varargin)
% PF0 M-file for pF0.fig
%      PF0, by itself, creates a new PF0 or raises the existing
%      singleton*.
%
%      H = PF0 returns the handle to a new PF0 or the handle to
%      the existing singleton*.
%
%      PF0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PF0.M with the given input arguments.
%
%      PF0('Property','Value',...) creates a new PF0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before processF0_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pF0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pF0

% Last Modified by GUIDE v2.5 12-May-2010 16:00:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pF0_OpeningFcn, ...
                   'gui_OutputFcn',  @pF0_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pF0 is made visible.
function pF0_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pF0 (see VARARGIN)

%Folders for the "other" processF0
rmpath('F:\neurostuff\2phAnalysis_pep\AnalysisCode')
rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI')
rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI\New')
rmpath('F:\neurostuff\2phAnalysis_pep\cbpep')
rmpath('F:\neurostuff\2phAnalysis_pep\pep2matlab_new')
rmpath('F:\neurostuff\2phAnalysis_pep\pepanalysis')

%Folders for the "other" processF0
rmpath('F:\neurostuff\2phAnalysis\AnalysisCode')
rmpath('F:\neurostuff\2phAnalysis\2ph_Processing')
rmpath('F:\neurostuff\2phAnalysis\2pAnGUI')
rmpath('F:\neurostuff\2phAnalysis\2pAnGUI\general')

%Folders for the "other" pF0
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\AnalysisCode')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\2pAnGUI')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\2pAnGUI\New')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\cbpep')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\pep2matlab_new')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\pepanalysis')
rmpath('Z:\Ian_N\Beta\2phAnalysis_pep\AnalysisCode\ContrastResp')

%Folders for this pF0
path('Z:\Ian_N\Beta\2phAnalysis\AnalysisCode',path)
path('Z:\Ian_N\Beta\2phAnalysis\2ph_Processing',path)
path('Z:\Ian_N\Beta\2phAnalysis\2pAnGUI',path)
path('Z:\Ian_N\Beta\2phAnalysis\2pAnGUI\general',path)
path('Z:\Ian_N\Beta\2phAnalysis\AnalysisCode\ContrastResp',path)

global processF0_handles

processF0_handles = handles;

% Choose default command line output for pF0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pF0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


clear f0m f0m_var Tens Tens_var repdom funcmap bcond bsflag bwCell1 symbolInfo 

% --- Outputs from this function are returned to the command line.
function varargout = pF0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
global f0m f0m_var Tens Tens_var AUE bsflag Flim repDom
% JHM 11/25/09 - f0m1 and f0m2 are for each channel a cell array of images where each
% image is the mean response for a condition within the defined time
% window. f0m1/2_var are variance for same time window across trial and
% time. Tens1/2 is a cell array, number of cells equals number of
% conditions, each cell contains an image stack where each slice is the
% average response at that time slice across trials of the same condition.
% Tens1/2_var is the same structure, but variance instead of mean. pepANA
% is the analyzer file information from the set directory. bsflag is the
% baseline subtraction flag (1 or 0). Flim is the delimeter for the time
% window. repDom vector that defines which trials to analyze (e.g. [1 3 5 7] or [1:2:end] to only analyze odd trials)

bsflag = 0;

t0 = cputime;

varflag = get(handles.trialVarianceFlag,'Value');
bsflag = get(handles.basesub,'Value');
Flim = str2double(get(handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(handles.epistop,'String')); %Frame stop in ms (to average)
b = str2double(get(handles.bstart,'String')); %in msec as well
b(2) = str2double(get(handles.bstop,'String'));

shiftFlag = get(handles.shiftFlag,'Value');  %lateral movement correction
lumFlag = get(handles.lumFlag,'Value');  %normalize luminance of each frame

repdum = get(handles.repDom,'string');  
if strcmp(repdum,'All')
    repDom = 1:getnorepeats(1);
else
    eval(['repDom = ' '[' repdum '];'])
end

setacqinfo(1)  %need to set to first trial so that it doesn't use a blank trial when doing CondTensor

set(handles.status,'string','Processing...'), drawnow

[Tens Tens_var] = CondTensor2(b,shiftFlag,lumFlag,varflag);  %%Compute entire space time block for each condition

f0m = CondF0(Tens,Flim);   %%%Compute mean over time interval [Flim(1) Flim(2)]%%%  
[f0m_var] = CondF0(Tens_var,Flim);

set(handles.status,'string','Done'), drawnow
%sound(.6*sin(2*pi*400/(1200)*(0:400)),1200)  %Signal done

t1 = cputime-t0;

set(handles.time,'string',num2str(t1))

set(handles.loaded,'string',AUE)

set(handles.plot,'enable','on')

set(handles.save,'enable','on')

set(handles.makeMask,'enable','on')



% --- Executes on selection change in func.
function func_Callback(hObject, eventdata, handles)
% hObject    handle to func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns func contents as cell array
%        contents{get(hObject,'Value')} returns selected item from func


% --- Executes during object creation, after setting all properties.
function func_CreateFcn(hObject, eventdata, handles)
% hObject    handle to func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
global f0m
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
    if isfield(S,'f1m')
        warndlg('This is processed data from an F1 experiment.  Try again.','!!!') 
    else
    f0m = S.f0cell{1};    %f0m is a cell array with images from each condition

    set(handles.plot,'enable','on')
    
    set(handles.loaded,'string',filename(1:length(filename)-4))
    
    end
end

function setimagedir_Callback(hObject, eventdata, handles)
% hObject    handle to setimagedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setimagedir as text
%        str2double(get(hObject,'String')) returns contents of setimagedir as a double


function loadexp_Callback(hObject, eventdata, handles)
% hObject    handle to loadexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadexp as text
%        str2double(get(hObject,'String')) returns contents of loadexp as a double


% --- Executes during object creation, after setting all properties.
function loadexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setdirs.
function setdirs_Callback(hObject, eventdata, handles)
% hObject    handle to setdirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ACQinfo twophDATADIR AUE bw Analyzer

anim = get(handles.loadana,'String');
expt = get(handles.loadexp,'String');

dir = get(handles.analyzedir,'String'); %partial path for analyzer file
setAnalyzerDirectory([dir anim '\'])

loadAnalyzer(expt)

fno = str2double(get(handles.frameno,'String'));    %Get frame number
tno = str2double(get(handles.trialno,'String'));    %Get trial number

dir = get(handles.datadir,'String'); %partial path for .tiff files 

twophDATADIR = [dir anim '\' expt '\'];  %Make path for .tiff files
AUE = [anim '_' expt]; %loaded unit experiment 'u000_000'

setacqinfo(tno)  %Set the global ACQinfo structure (contains trial info)

[Im] = Load2phImage(fno,[1 1 0],tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)       %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])

conds = getnoconditions;
reps = getnorepeats(1);
set(handles.nocond,'string',num2str(conds))
set(handles.norep,'string',num2str(reps))
set(handles.dirstatus,'string','Loaded')

set(handles.setROI,'enable','on')
set(handles.process,'enable','on')

if isfield(ACQinfo,'stimPredelay')  %this wasn't saved in the earliest experiments (aa0 to aa2)
    set(handles.predelay,'string',['predelay=' num2str(getParamVal('predelay'))])
    set(handles.postdelay,'string',['postdelay=' num2str(getParamVal('postdelay'))])
    set(handles.trialtime,'string',['trialtime=' num2str(getParamVal('stim_time'))])
end
        
Nsym = length(Analyzer.loops.conds{1}.symbol);

set(handles.primSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
set(handles.primSymbol,'value',1)
set(handles.secSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
set(handles.secSymbol,'value',2)
set(handles.tertSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
set(handles.tertSymbol,'value',3)

switch Nsym   
    case 1        
        set(handles.secSymbol,'enable','off')
        set(handles.tertSymbol,'enable','off')
        set(handles.secCollapse,'enable','off')
        set(handles.tertCollapse,'enable','off')        
    case 2        
        set(handles.secSymbol,'enable','on')
        set(handles.tertSymbol,'enable','off')
        set(handles.secCollapse,'enable','on')
        set(handles.tertCollapse,'enable','off')                
    case 3        
        set(handles.secSymbol,'enable','on')
        set(handles.tertSymbol,'enable','on')
        set(handles.secCollapse,'enable','on')
        set(handles.tertCollapse,'enable','on')
end

if ~isempty(bw)
    if length(bw(:,1)) ~= length(Im{1}(:,1)) || length(bw(1,:)) ~= length(Im{1}(1,:))
        bw = ones(size(Im{1}));
    end
else
    bw = ones(size(Im{1}));
end

function datadir_Callback(hObject, eventdata, handles)
% hObject    handle to datadir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datadir as text
%        str2double(get(hObject,'String')) returns contents of datadir as a double


% --- Executes during object creation, after setting all properties.
function datadir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datadir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epistart_Callback(hObject, eventdata, handles)
% hObject    handle to epistart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epistart as text
%        str2double(get(hObject,'String')) returns contents of epistart as a double


% --- Executes during object creation, after setting all properties.
function epistart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epistart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epistop_Callback(hObject, eventdata, handles)
% hObject    handle to epistop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epistop as text
%        str2double(get(hObject,'String')) returns contents of epistop as a double


% --- Executes during object creation, after setting all properties.
function epistop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epistop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HPBW_Callback(hObject, eventdata, handles)
% hObject    handle to HPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HPBW as text
%        str2double(get(hObject,'String')) returns contents of HPBW as a double


% --- Executes during object creation, after setting all properties.
function HPBW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LPBW_Callback(hObject, eventdata, handles)
% hObject    handle to LPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LPBW as text
%        str2double(get(hObject,'String')) returns contents of LPBW as a double


% --- Executes during object creation, after setting all properties.
function LPBW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bstart_Callback(hObject, eventdata, handles)
% hObject    handle to bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bstart as text
%        str2double(get(hObject,'String')) returns contents of bstart as a double


% --- Executes during object creation, after setting all properties.
function bstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bstop_Callback(hObject, eventdata, handles)
% hObject    handle to bstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bstop as text
%        str2double(get(hObject,'String')) returns contents of bstop as a double


% --- Executes during object creation, after setting all properties.
function bstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in basesub.
function basesub_Callback(hObject, eventdata, handles)
% hObject    handle to basesub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of basesub

% --- Executes on button press in tempfilt.
function tempfilt_Callback(hObject, eventdata, handles)
% hObject    handle to tempfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tempfilt


function Hwidth_Callback(hObject, eventdata, handles)
% hObject    handle to Hwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hwidth as text
%        str2double(get(hObject,'String')) returns contents of Hwidth as a double


% --- Executes during object creation, after setting all properties.
function Hwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in HPWind.
function HPWind_Callback(hObject, eventdata, handles)
% hObject    handle to HPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns HPWind contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HPWind


% --- Executes during object creation, after setting all properties.
function HPWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lwidth_Callback(hObject, eventdata, handles)
% hObject    handle to Lwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lwidth as text
%        str2double(get(hObject,'String')) returns contents of Lwidth as a double


% --- Executes during object creation, after setting all properties.
function Lwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LPWind.
function LPWind_Callback(hObject, eventdata, handles)
% hObject    handle to LPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LPWind contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LPWind


% --- Executes during object creation, after setting all properties.
function LPWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HPflag.
function HPflag_Callback(hObject, eventdata, handles)
% hObject    handle to HPflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HPflag


% --- Executes on button press in LPflag.
function LPflag_Callback(hObject, eventdata, handles)
% hObject    handle to LPflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LPflag


% --- Executes on button press in setROI.
function setROI_Callback(hObject, eventdata, handles)
global bw f0m
% hObject    handle to setROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fno = str2double(get(handles.frameno,'String'));    %Get frame number
tno = str2double(get(handles.trialno,'String'));    %Get frame number

[Im] = Load2phImage(fno,[1 1 0],tno);

figure,imagesc(Im{1}), colormap gray        

bw = roipoly;
close

if ~isempty(f0m)
    set(handles.plot,'enable','on')
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
global bw f0m funcmap bcond ACQinfo bwCell1 TCWin symbolInfo Analyzer
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(bw) || length(bw(:,1)) ~= ACQinfo.linesPerFrame || length(bw(1,:)) ~= ACQinfo.pixelsPerLine
    bw = ones(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
end

TCWin = str2double(get(handles.TCWin,'string'));

checkvect = [get(handles.F0im,'value') get(handles.funcim,'value')];

togstateHP = get(handles.HPflag,'Value');
togstateLP = get(handles.LPflag,'Value');

sizeIM = size(f0m{1});
if togstateHP == 1
    Hwidth = str2double(get(handles.Hwidth,'string'));
    ind = get(handles.HPWind,'value');

    switch ind
        case 1
            H = -fspecial('gaussian',sizeIM,Hwidth);
            H(round(sizeIM/2),round(sizeIM/2)) = 1+H(round(sizeIM/2),round(sizeIM/2));
        case 2
            H = zeros(sizeIM);
            Hd = hann(Hwidth)*hann(Hwidth)';
            Hd = -Hd./sum(Hd(:));
            Hd(round(Hwidth/2),round(Hwidth/2)) = 1+Hd(round(Hwidth/2),round(Hwidth/2));
            H(1:Hwidth,1:Hwidth) = Hd;
        case 3
            H = zeros(sizeIM);
            Hd = -fspecial('disk',round(Hwidth/2));
            Hsize = length(Hd(1,:));  %~=Hwidth
            Hd(round(Hsize/2),round(Hsize/2)) = 1+Hd(round(Hsize/2),round(Hsize/2));
            H(1:Hsize,1:Hsize) = Hd;
    end
    if togstateLP == 0
        hh = ifft2(abs(fft2(H)));   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(handles.Lwidth,'string'));
    ind = get(handles.LPWind,'value');

    switch ind
        case 1
            L = fspecial('gaussian',sizeIM,Lwidth);
        case 2
            L = zeros(sizeIM);
            Ld = hann(Lwidth)*hann(Lwidth)';
            Ld = Ld./sum(Ld(:));
            L(1:Lwidth,1:Lwidth) = Ld;
        case 3
            L = zeros(sizeIM);
            Ld = fspecial('disk',round(Lwidth/2));
            Lsize = length(Ld(1,:));
            L(1:Lsize,1:Lsize) = Ld;
    end
    if togstateHP == 0
        hh = ifft2(abs(fft2(L)));   %Eliminate phase information
    else
        hh = ifft2(abs(fft2(L).*fft2(H)));   %Take mag because phase gives a slight shift.
    end
end

if ~or(togstateLP,togstateHP)
    hh = [];
end

%%...Done making filter

anatomyflag = get(handles.anatomyFlag,'value');
cellmaskflag = get(handles.cellmaskflag,'value');

%%%%%%Put all the symbol information into global structure\
symbolInfo = struct;
Fsymbol = get(handles.primSymbol,'string'); %primary parameter symbol in looper to analyze
symbolInfo.ID(1) = get(handles.primSymbol,'value');  %The index with respect to the looper
symbolInfo.str{1} = Fsymbol{symbolInfo.ID(1)};  %Selected string
symbolInfo.domType = get(handles.domType,'value');  %Type of domain for primary symbol... .e.g. circular 'Axis'

Fsymbol = get(handles.secSymbol,'string'); %secondary symbol
symbolInfo.ID(2) = get(handles.secSymbol,'value');
symbolInfo.str{2} = Fsymbol{symbolInfo.ID(2)};
symbolInfo.Collapse(1) = get(handles.secCollapse,'value');  %Describes how to collapse across secondary loop domains

Fsymbol = get(handles.tertSymbol,'string'); %tertiary symbol
symbolInfo.ID(3) = get(handles.tertSymbol,'value');
symbolInfo.str{3} = Fsymbol{symbolInfo.ID(3)};
symbolInfo.Collapse(2) = get(handles.tertCollapse,'value');  %Describes how to collapse across tertiary loop domains
%%%%%%%%%%%%

%%Filter raw F0 images with hh and create the functional maps...
if checkvect(2) == 1        %if "Functional Images" is checked
    
%     if cellmaskflag && ~isempty(bwCell1) 
%         [kernPop popBlank CoM pardom] = GetROIKernels(f0m,bwCell1);
%     end
    
    switch symbolInfo.domType

        case 1        %domain type is 'axis'
            funcmap = GprocessAxis(f0m,hh);  %output is a vector image

        case 2        %domain type is 'direction'
            funcmap = GprocessDir(f0m,hh);  %output is a vector image

        case 3        %domain type is 'log'
            funcmap = GprocessLog(f0m,bw,hh);  %output is complex image: real(funcmap) = mag; imag(funcmap) = pref

        case 4        %domain type is '2D'
            [angx magx angy magy] = Gprocessret(f0m,hh);
    end

%     if funcflag == 6        %%functionality is color
%         [funcmap colorSel] = GprocessColor(f0m,hh);
%     end
%     
%     if funcflag == 7        %%functionality is color & form
%         [funcmap colorSel] = GprocessColorForm(f0m,hh);
%     end
end


%If the radio button is clicked AND a mask exists,
%then average pixels within all the ROIs... "cellurize"
if cellmaskflag && ~isempty(bwCell1)  
    funcmap = Cellurize(funcmap,bwCell1);    
    bwCellPlot = bwCell1;
else
    bwCellPlot = ones(size(funcmap));
end

%Create plots

switch symbolInfo.domType

    case 1    %"Axis"

        if checkvect(1) == 1  %mean image for each condition
            plotF0(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            ang = angle(funcmap);
            mag = abs(funcmap);
            ang = (ang+pi*(1-sign(ang)))/2*180/pi;  %Put domain as [0 180].
            figure            
            Gplotaxismap(mag.*bwCellPlot.*bw,ang,anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);

            %GplotorimapDots(ang,celllocs1)   %To use this, other things must change
        end

    case 2    %"direction"
        
        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            ang = angle(funcmap);
            mag = abs(funcmap);
            ang = (ang+pi*(1-sign(ang)))*180/pi;  %make it 0 to 360.
            figure
            Gplotdirmap(mag.*bwCellPlot.*bw,ang,anatomyflag), title('Direction','FontWeight','bold','FontSize',15);
                
            %GplotorimapDots(ang,celllocs1)   %To use this other things
            %must change
        end

    case 3    %"log domain"
        
        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            
            Int = real(funcmap);
            pref = imag(funcmap);
            figure
            Gplotlogmap(Int.*bwCellPlot.*bw,pref,anatomyflag), title(symbolInfo.str{1}(find(symbolInfo.str{1}~='_')),'FontWeight','bold','FontSize',15);

        end

    case 4    %Retinotopy
        
        if checkvect(1) == 1
            
            mima = [min(f0m{1}(:)) max(f0m{1}(:))];
            N = length(f0m);
            k = 1;
            figure
            for i = 1:N
                if i ~= bcond+1
                    subplot(1,N-length(bcond),k)
                    imagesc(f0m{i},mima)
                    title(['Condition ' num2str(k-1)])

                    k = k+1;
                end
            end
            colormap gray
        end
        
        if checkvect(2) == 1  %functional maps

            screenDist = Analyzer.M.screenDist;
            screenResX = Analyzer.M.xpixels/Analyzer.M.screenXcm;  %pix/cm
            screenResY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

            [xpos ypos xsize ysize] = getPosSize;
            xsize_cm = (max(xpos)-min(xpos)+min(xsize))/screenResX;  %cm stimulus width
            %xsize_deg = 2*atan2(xsize_cm/2,screenDist)*180/pi;  %convert to deg
            xsize_deg = 360*xsize_cm/(2*pi*screenDist);
            
            ysize_cm = (max(ypos)-min(ypos)+min(ysize))/screenResY;  %cm stimulus width
            %ysize_deg = 2*atan2(ysize_cm/2,screenDist)*180/pi;  %convert to deg
            ysize_deg = 360*ysize_cm/(2*pi*screenDist);

            figure('Name','RETINOTOPIC MAPS','NumberTitle','off')
            subplot(2,1,1), imagesc(angx,'AlphaData',magx.*bw,[0 xsize_deg]); colorbar, set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]), ylabel('x position (deg)')
            subplot(2,1,2), imagesc(angy,'AlphaData',magy.*bw,[0 ysize_deg]); colorbar, set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]), ylabel('y position (deg)')
            plotpixel_cb

            figure('Name','SCREEN COVERAGE','NumberTitle','off')
            Gplotretcoverage(angx,magx,angy,magy)

        end

%     case 6    %"color"
%         
%         if checkvect(1) == 1
%             plotF0(f0m,bw,hh)
%         end
%         if checkvect(2) == 1  %functional maps
%             %funcmap = Cellurize(funcmap,bwCell1);
%             ang = angle(funcmap);
%             mag = abs(funcmap);
%             ang = (ang+pi*(1-sign(ang)))/2*180/pi;  %Put in orientation domain and convert to degrees.
%             figure
%             Gplotcolormap(bw.*mag,ang), title('Color Axis','FontWeight','bold','FontSize',15);
% 
%             figure
%             subplot(1,2,1)
%             imagesc(colorSel), colormap gray, colorbar
%             colorbar('YTick',[-1 0 1],'YTickLabel',{'-1 Lum','0 Both','+1 Color'})
%             title('Color to Lum Response','FontWeight','bold','FontSize',12);
% 
%         end

end


% --- Executes on button press in F0im.
function F0im_Callback(hObject, eventdata, handles)
% hObject    handle to F0im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F0im



% --- Executes on button press in funcim.
function funcim_Callback(hObject, eventdata, handles)
% hObject    handle to funcim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of funcim



function frameno_Callback(hObject, eventdata, handles)
% hObject    handle to frameno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameno as text
%        str2double(get(hObject,'String')) returns contents of frameno as a double

tno = str2double(get(handles.trialno,'String'));    %Get trial number
fno = str2double(get(handles.frameno,'String'));    %Get frame number
[Im] = Load2phImage(fno,[1 1 0],tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])

% --- Executes during object creation, after setting all properties.
function frameno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
global f0m
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UE = get(handles.loadexp,'string');
path = 'c:\Processed Data\';
filename = strcat(path,AUE);
f0cell{1} = f0m;   %Channel 1
uisave('f0cell',filename)


function trialno_Callback(hObject, eventdata, handles)
% hObject    handle to trialno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trialno as text
%        str2double(get(hObject,'String')) returns contents of trialno as a double

tno = str2double(get(handles.trialno,'String'));    %Get trial number
fno = str2double(get(handles.frameno,'String'));    %Get frame number
[Im] = Load2phImage(fno,[1 1 0],tno);

setacqinfo(tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot framesetacqinfo(tno);  %Set global ACQinfo for this trial
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])

% --- Executes during object creation, after setting all properties.
function trialno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function analyzedir_Callback(hObject, eventdata, handles)
% hObject    handle to analyzedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analyzedir as text
%        str2double(get(hObject,'String')) returns contents of analyzedir as a double


% --- Executes during object creation, after setting all properties.
function analyzedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyzedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loadana_Callback(hObject, eventdata, handles)
% hObject    handle to loadana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadana as text
%        str2double(get(hObject,'String')) returns contents of loadana as a double


% --- Executes during object creation, after setting all properties.
function loadana_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in lumFlag.
function lumFlag_Callback(hObject, eventdata, handles)
% hObject    handle to lumFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lumFlag


% --- Executes on button press in shiftFlag.
function shiftFlag_Callback(hObject, eventdata, handles)
% hObject    handle to shiftFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shiftFlag




% --- Executes on button press in cellmaskflag.
function cellmaskflag_Callback(hObject, eventdata, handles)
% hObject    handle to cellmaskflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cellmaskflag




% --- Executes on button press in anatomyFlag.
function anatomyFlag_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anatomyFlag




% --- Executes on button press in cellmaskflag.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to cellmaskflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cellmaskflag


% --- Executes on button press in anatomyFlag.
function radiobutton15_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anatomyFlag



function maskSize_Callback(hObject, eventdata, handles)
% hObject    handle to maskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskSize as text
%        str2double(get(hObject,'String')) returns contents of maskSize as a double


% --- Executes during object creation, after setting all properties.
function maskSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskThresh_Callback(hObject, eventdata, handles)
% hObject    handle to maskThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskThresh as text
%        str2double(get(hObject,'String')) returns contents of maskThresh as a double


% --- Executes during object creation, after setting all properties.
function maskThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskMorph_Callback(hObject, eventdata, handles)
% hObject    handle to maskMorph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskMorph as text
%        str2double(get(hObject,'String')) returns contents of maskMorph as a double


% --- Executes during object creation, after setting all properties.
function maskMorph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskMorph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in trialVarianceFlag.
function trialVarianceFlag_Callback(hObject, eventdata, handles)
% hObject    handle to trialVarianceFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trialVarianceFlag





function repDom_Callback(hObject, eventdata, handles)
% hObject    handle to repDom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repDom as text
%        str2double(get(hObject,'String')) returns contents of repDom as a double


% --- Executes during object creation, after setting all properties.
function repDom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repDom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in openPOPgui.
function openPOPgui_Callback(hObject, eventdata, handles)
% hObject    handle to openPOPgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PopAnalysis


% --- Executes on button press in makeMask.
function makeMask_Callback(hObject, eventdata, handles)
% hObject    handle to makeMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global bwCell1 bwCell2

[bwCell1 bwCell2] = MakeCellMask(str2num(get(handles.maskSize,'string')),str2num(get(handles.maskThresh,'string')),str2num(get(handles.maskMorph,'string')),str2num(get(handles.minCellArea,'string')));

set(handles.openPOPgui,'enable','on')



function searchRange_Callback(hObject, eventdata, handles)
% hObject    handle to searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of searchRange as text
%        str2double(get(hObject,'String')) returns contents of searchRange as a double


% --- Executes during object creation, after setting all properties.
function searchRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function minCellArea_Callback(hObject, eventdata, handles)
% hObject    handle to minCellArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minCellArea as text
%        str2double(get(hObject,'String')) returns contents of minCellArea as a double


% --- Executes during object creation, after setting all properties.
function minCellArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minCellArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function TCWin_Callback(hObject, eventdata, handles)
% hObject    handle to TCWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TCWin as text
%        str2double(get(hObject,'String')) returns contents of TCWin as a double


% --- Executes during object creation, after setting all properties.
function TCWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TCWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function primSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to primSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of primSymbol as text
%        str2double(get(hObject,'String')) returns contents of primSymbol as a double


% --- Executes during object creation, after setting all properties.
function primSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to primSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in domType.
function domType_Callback(hObject, eventdata, handles)
% hObject    handle to domType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns domType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from domType


% --- Executes during object creation, after setting all properties.
function domType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to domType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in secCollapse.
function secCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to secCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns secCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secCollapse


% --- Executes during object creation, after setting all properties.
function secCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in secSymbol.
function secSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to secSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns secSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secSymbol


% --- Executes during object creation, after setting all properties.
function secSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tertCollapse.
function tertCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to tertCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tertCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tertCollapse


% --- Executes during object creation, after setting all properties.
function tertCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tertCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tertSymbol.
function tertSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to tertSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tertSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tertSymbol


% --- Executes during object creation, after setting all properties.
function tertSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tertSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


