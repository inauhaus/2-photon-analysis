function varargout = pF0(varargin)

%Ian Nauhaus

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

% Last Modified by GUIDE v2.5 06-May-2020 09:54:38

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
% rmpath('F:\neurostuff\2phAnalysis_pep\AnalysisCode')
% rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI')
% rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI\New')
% rmpath('F:\neurostuff\2phAnalysis_pep\cbpep')
% rmpath('F:\neurostuff\2phAnalysis_pep\pep2matlab_new')
% rmpath('F:\neurostuff\2phAnalysis_pep\pepanalysis')
% 
% 
% %Folders for the "other" processF0
% rmpath('F:\neurostuff\2phAnalysis\AnalysisCode')
% rmpath('F:\neurostuff\2phAnalysis\2ph_Processing')
% rmpath('F:\neurostuff\2phAnalysis\AnalysisCode\DynamicProcess')
% rmpath('F:\neurostuff\2phAnalysis\2pAnGUI')
% rmpath('F:\neurostuff\2phAnalysis\2pAnGUI\general')
% 
% %Folders for the "other" pF0
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\2pAnGUI')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\2pAnGUI\New')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\cbpep')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\pep2matlab_new')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\pepanalysis')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\ContrastResp')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\Masking')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\DynamicProcess')

p = 'C:\Users\Ian Nauhaus\Documents\Matlab\2pScanboxAnalysis\';
path([p '\AnalysisCode'],path)
path([p '\KalRet'],path)
path([p '\2ph_Processing'],path)
path([p '\2pAnGUI'],path)
path([p '\2pAnGUI\general'],path)
path([p '\AnalysisCode\ContrastResp'],path)
path([p '\AnalysisCode\KalRet'],path)
path([p '\AnalysisCode\DynamicProcess'],path)
path([p '\AnalysisCode\DynamicProcess\RevCorr_GUI'],path)
path([p '\AnalysisCode\Masking'],path)
path([p '\AnalysisCode\Color'],path)
path([p '\AnalysisCode\Tensor_processing'],path)
path([p '\offlineMovementCorrection'],path)
path([p '\monkeyMaps'],path)
path([p '\monkeyMaps'],path)
path([p '\Model Fitting'],path)
path([p '\mouse_color'],path)

path('C:\sbx_files\sbx',path)



global G_handles

G_handles = handles;

% Choose default command line output for pF0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pF0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


clear global f0m f0m_var Tens Tens_var repdom funcmap bcond bsflag maskS symbolInfo 

% --- Outputs from this function are returned to the command line.
function varargout = pF0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in computeF0.
function computeF0_Callback(hObject, eventdata, handles)

global Analyzer

processButton(1);

% --- Executes on button press in computeTraces.
function computeTraces_Callback(hObject, eventdata, handles)
% hObject    handle to computeTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

if isfield(maskS,'bwCell')
    processButton(0);
else
    'Need to load mask first'
end

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


% --- Executes on button press in loadF0.
function loadF0_Callback(hObject, eventdata, handles)

% hObject    handle to loadF0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global f0m Tens 

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
    if isfield(S,'f0m')
        f0m = S.f0m;    %f0m is a cell array with images from each condition
        Tens = S.Tens;
    end

    set(handles.plot,'enable','on')
    
    set(handles.loaded,'string',filename(1:length(filename)-4))

end

% --- Executes on button press in loadTraces.
function loadTraces_Callback(hObject, eventdata, handles)
% hObject    handle to loadTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cellS maskS

[filename, pathname] = uigetfile({'*.mat;*.signals'}, 'Pick a file');

if filename ~= 0
    id = find(filename == '.');
    ext = filename(id+1:end);
    
    if strcmp(ext,'mat')
        
        S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
        
        cellS = S.cellS;
        set(handles.plot,'enable','on')
        
        set(handles.loaded,'string',filename(1:length(filename)-4))
        
    elseif strcmp(ext,'signals') %Via scanbox cell extractor
        
        %Load the matrix of signals
        S = load(strcat(pathname,filename),'-mat');  %Returns the contents in the .mat under the structure S
        set(handles.plot,'enable','on')        
        set(handles.loaded,'string',filename(1:length(filename)-4))        
        
        %Load corresponding cell mask from Scanbox
        a = questdlg('Load the corresponding cell mask created in scanbox?');
        if strcmp('Yes',a)            
            [filename, pathname] = uigetfile('*.segment', 'Pick a .mat or .segment file');
            loadMaskInfo_SBX(pathname,filename)
            figure(40),
            %imagesc(maskS.im{1}), colormap gray
            %hold on
            contour(maskS.bwCell{1},[.5 .5],'r')           
        end
        
        %Reorganize timecourses        
        signew = zeros([size(S.sig,1) size(S.sig,2)+1]);
        idim = bwlabel(sign(maskS.bwCell{1}));
        cellID = unique(idim);
        %This leaves the initial entry (i.e. neuropil) of signew = 0. 
        for i = 2:length(cellID)           
            idx = find(idim == cellID(i));      
            sbxid = maskS.bwCell_sbxID{1}(idx(1));
            signew(:,i) = S.sig(:,sbxid);
        end       
        
        %Create cellS structure
        parseToCellS(signew) %Create trial blocks: cellS.cellMat
        getCellStats
                
    end
    
end

function setimagedir_Callback(hObject, eventdata, handles)
% hObject    handle to setimagedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setimagedir as text
%        str2double(get(hObject,'String')) returns contents of setimagedir as a double


% --- Executes on button press in setdirs.
function setdirs_Callback(hObject, eventdata, handles)
% hObject    handle to setdirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ACQinfo Analyzer G_handles cellS

Gsetdirectories %Load analyzer file and scanbox info

fno = 1; tno = 1;
set(handles.frameno,'string',num2str(fno))
set(handles.trialno,'string',num2str(tno))
[Im] = Load2phImage(fno,tno);

axes(G_handles.rimage1);     %Make rimage current figure
cla
imageRGB(double(Im{1}),2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(double(Im{2}),1)       %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])

setGUIlabels

%handle the slider
nF = median(diff(ACQinfo.SBInfo.frame(2:end-1))); %frames/trial
nT = getnotrials;
set(G_handles.frameslide,'value',(fno-1)/(nF-1))
set(G_handles.trialslide,'value',(tno-1)/(nT-1))

stepsize = 1/nF;
set(G_handles.frameslide,'SliderStep',[stepsize 4*stepsize])
stepsize = 1/nT;
set(G_handles.trialslide,'SliderStep',[stepsize 4*stepsize])

set(G_handles.makeTemplate,'enable','on')

fixSBsyncs

anatomyTemp_Callback  %Creates maskS.anatomyTemplate;

cellS = struct; %clear it

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

[Im] = Load2phImage(fno,tno);

figure,imagesc(Im{1}), colormap gray        

bw = roipoly;
close

if ~isempty(f0m)
    set(handles.plot,'enable','on')
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
global bw f0m funcmap bcond ACQinfo maskS TCWin symbolInfo Analyzer
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(bw) || length(bw(:,1)) ~= ACQinfo.linesPerFrame || length(bw(1,:)) ~= ACQinfo.pixelsPerLine
    bw = ones(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
end

TCWin = str2double(get(handles.TCWin,'string'));

checkvect = [get(handles.F0im,'value') get(handles.funcim,'value')];

hh = makeMapFilter;

anatomyflag = get(handles.anatomyFlag,'value');

setsymbolstruct %Put all the symbol information into global structure

%%Filter raw F0 images with hh and create the functional maps...
if checkvect(2) == 1        %if "Functional Images" is checked
    
    switch symbolInfo.domType

        case 1        %domain type is 'axis'
            funcmap = GprocessAxis(f0m,hh);  %output is a vector image

        case 2        %domain type is 'direction'
            funcmap = GprocessDir(f0m,hh);  %output is a vector image

        case 3        %domain type is 'log'
            funcmap = GprocessLog(f0m,bw,hh);  %output is complex image: real(funcmap) = mag; imag(funcmap) = pref

        case 4        %domain type is '2D'
            [angx magx angy magy] = Gprocessret(f0m,hh);
            
        case 5        %domain type is 'binary'
            funcmap = GprocessBinary(f0m,hh);
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
if get(handles.cellmaskplotflag,'value') && ~isempty(maskS.bwCell{1})  

    getCellStats;  %reset the the mean based on time window
    getNeuronMask;
    funcmap = Cellurize(funcmap,maskS.neuronmask);
    bwCellPlot = maskS.neuronmask;
    %funcmap = Cellurize(funcmap,maskS.bwCell{1});    
    %bwCellPlot = maskS.bwCell{1};
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
            if get(handles.cellmaskplotflag,'value')
                Gplotaxismap_cellmask(mag.*bwCellPlot.*bw,ang,anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);
                %Gplotaxismap_cellmask_Gfit(mag.*bwCellPlot.*bw,ang,anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);
            else
                Gplotaxismap(mag.*bwCellPlot.*bw,ang,anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);
            end

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
        
    case 5    %"binary domain" (e.g. ocular dominance)

        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        
        if checkvect(2) == 1  %functional maps
            mi = prctile(funcmap(:),2); ma = prctile(funcmap(:),98);
            figure, imagesc(funcmap,[mi ma]), colormap gray
            title(symbolInfo.str{1}(find(symbolInfo.str{1}~='_')),'FontWeight','bold','FontSize',15);
            axis square
        end

%     case 6    %"color"
%         
%         if checkvect(1) == 1
%             plotF0(f0m,bw,hh)
%         end
%         if checkvect(2) == 1  %functional maps
%             %funcmap = Cellurize(funcmap,maskS.bwCell{1});
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



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)

% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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


% --- Executes on button press in lumFlag.
function lumFlag_Callback(hObject, eventdata, handles)
% hObject    handle to lumFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lumFlag


% --- Executes on button press in slowXcorr.
function slowXcorr_Callback(hObject, eventdata, handles)
% hObject    handle to slowXcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slowXcorr



% --- Executes on button press in anatomyFlag.
function anatomyFlag_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anatomyFlag



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




% --- Executes on button press in EbarFlag.
function EbarFlag_Callback(hObject, eventdata, handles)
% hObject    handle to EbarFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EbarFlag


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


% --- Executes on button press in makeTemplate.
function makeTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to makeTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS f0m ACQinfo

%'f0forMask' uses the entire experiment to compute the template for the mask.
%N.B. computef0.raw must have already been run to compute the 'f0m' images...
%and w/o baseline subtraction
%%%%%%

Tdelim = str2num(get(handles.maskTemplateDelim,'string'));
fp = ACQinfo.linesPerFrame*ACQinfo.msPerLine/1000;
predelayFrames = round(getparam('predelay')/fp);

im = 0;
global twophDATADIR
fr = ACQinfo.SBInfo.resfreq/ACQinfo.SBInfo.sz(1); %frames/sec
for i = 1:length(Tdelim)-1
    Twin = [Tdelim(i) Tdelim(i+1)];
    Fwin = round(Twin*fr)+1;
    dumdum = double(squeeze(sbxread(twophDATADIR,Fwin(1),1+Fwin(2)-Fwin(1))));
    dum{1} = dumdum(:,ACQinfo.unblanked,:);
    
    
    switch get(handles.templateType,'value')
        case 1
            dum{1} = cleanTensor4(dum{1},Fwin); %This will downsample it in space and time
            im = im + getMaxProj(dum{1},500)/(length(Tdelim)-1);            
        case 2
            dum{1} = cleanTensor4(dum{1},Fwin); %This will downsample it in space and time
            im = im + getSkewProj(dum{1},500)/(length(Tdelim)-1);
        case 3
            dum{1} = cleanTensor4(dum{1},Fwin); %This will downsample it in space and time
            im = im + mean(dum{1},3)/(length(Tdelim)-1);
        case 4
            dum{1} = cleanTensor4(dum{1},Fwin); %This will downsample it in space and time
            im = im + median(dum{1},3)/(length(Tdelim)-1);
        case 5
            
            %dum{1} = cleanTensor3(dum{1},Fwin);
            dum{1} = cleanTensor4(dum{1},Fwin,2); %This will downsample it in space and time.
                   
            %im = im + localXCorr2(dum{1},2)/(length(Tdelim)-1);
            %im = im + localXCorr(dum{1})/(length(Tdelim)-1);      
            im = im + localXCorr3(dum{1},2)/(length(Tdelim)-1);
            %im = im + localXCorr4(dum{1},3);
    end
    
    %im = LocalZ(im,100,1);
    maskS.im{1} = im;
        
end


maskS.bwCell{1} = zeros(length(ACQinfo.ydom),length(ACQinfo.xdom));
% if ~isempty('maskS')
%     if isfield(maskS,'bwCell');
%         if ~isempty(maskS.bwCell{1})
%             maskS.bwCell{1} = zeros(length(ACQinfo.ydom),length(ACQinfo.xdom));
%         end
%     end
% end

dim = size(maskS.im{1});
dimI = [ACQinfo.SBInfo.sz(1) length(ACQinfo.unblanked)];
maskS.im{1} = interp1(1:dim(1),maskS.im{1},linspace(1,dim(1),dimI(1)));
maskS.im{1} = interp1(1:dim(2),maskS.im{1}',linspace(1,dim(2),dimI(2)))';

maskS.bw = zeros(size(maskS.im{1}));


[xmicperpix ymicperpix] = getImResolution;
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
figure(40), 
maskS.im{1} = (maskS.im{1} - mean(maskS.im{1}(:)))/std(maskS.im{1}(:));
maskS.im{1}(find(maskS.im{1}<-2)) = -2;
maskS.im{1}(find(maskS.im{1}>6)) = 6;
imagesc(xdom,ydom,maskS.im{1}), colormap gray
axis image
hold off

%maskS.bwCell{1} = zeros(size(maskS.im{1}));

set(handles.saveMask,'enable','on')



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

global Analyzer G_handles

symid = get(G_handles.secSymbol,'value');

dom = getdomain(Analyzer.loops.conds{1}.symbol{symid});
domstr{1} = 'mean';
for i = 2:length(dom)+1
    domstr{i} = dom(i-1);
end
set(G_handles.secCollapse,'string',domstr)

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

global Analyzer G_handles

symid = get(G_handles.tertSymbol,'value');

dom = getdomain(Analyzer.loops.conds{1}.symbol{symid});
domstr{1} = 'mean';
for i = 2:length(dom)+1
    domstr{i} = dom(i-1);
end
set(G_handles.secCollapse,'string',domstr)

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




% --- Executes on slider movement.
function trialslide_Callback(hObject, eventdata, handles)
% hObject    handle to trialslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Nt = getnotrials;
tno = get(handles.trialslide,'value');
tno = round((Nt-1)*tno)+1;

set(handles.trialno,'string',num2str(tno));    %Get trial number
fno = str2double(get(handles.frameno,'string'));    %Get frame number
[Im] = Load2phImage(fno,tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(double(Im{1}),2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(double(Im{2}),1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])


% --- Executes during object creation, after setting all properties.
function trialslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function frameslide_Callback(hObject, eventdata, handles)
% hObject    handle to frameslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global ACQinfo

tno = str2double(get(handles.trialno,'string'));    %Get trial number

fno = get(handles.frameslide,'value');

frame1 = ACQinfo.SBInfo.frame((tno*2)-1); %first frame of trial is every other sync
Nf =  length(frame1:ACQinfo.SBInfo.frame((tno*2)));

fnoTrial = round((Nf-1)*fno)+1; %Within trial index
%fnoExp = fnoTrial + frame1; %Within experiment index.

set(handles.frameno,'string',num2str(fnoTrial)); %Show frame number in trial

[Im] = Load2phImage(fnoTrial,tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(double(Im{1}),2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(double(Im{2}),1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])


% --- Executes during object creation, after setting all properties.
function frameslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in fastXcorr.
function fastXcorr_Callback(hObject, eventdata, handles)
% hObject    handle to fastXcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fastXcorr


% --- Executes on selection change in motionModel.
function motionModel_Callback(hObject, eventdata, handles)
% hObject    handle to motionModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns motionModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from motionModel


% --- Executes during object creation, after setting all properties.
function motionModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to motionModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in saveMask.
function saveMask_Callback(hObject, eventdata, handles)
% hObject    handle to saveMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global maskS

path = 'c:\CellMasks\';
filename = path
uisave({'maskS'},filename)


% --- Executes on button press in loadMask.
function loadMask_Callback(hObject, eventdata, handles)
% hObject    handle to loadMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global maskS

anatdum = maskS.anatomyTemplate;  %save this, because loading the mask will save over it.

[filename, pathname] = uigetfile('*.mat;*.segment', 'Pick a .mat or .segment file');


if filename ~= 0
    
    ext = filename(find(filename == '.')+1:end);
    
    if strcmp(ext,'mat')
        
        loadMaskInfo(pathname,filename)
        figure(40),
        imagesc(maskS.im{1}), colormap gray
        hold on
        contour(maskS.bwCell{1},[.5 .5],'r')
        
        maskS.bwCell_sbxID{1} = []; %Make sure this gets cleared
        
    elseif strcmp(ext,'segment')
        
        loadMaskInfo_SBX(pathname,filename)
        figure(40),
        %imagesc(maskS.im{1}), colormap gray
        %hold on
        contour(maskS.bwCell{1},[.5 .5],'r')
        
        
    end
end

temp_ExptX = maskS.anatomyTemplate; %save this for alignment purposes below
maskS.anatomyTemplate = anatdum; %Put back the template for the loaded experiment (generated when hitting set dirs)


a = questdlg('Is this mask is from a different experiment and therefore be aligned to the loaded experiment directory?');
if strcmp('Yes',a)

    %Get shift values to align the loaded mask expt to the experiment in
    %the current directory. The loaded mask is aligned to temp_ExptX.  We
    %want it aligned to maskS.anatomyTemplate
    [mbest nbest] = getShiftVals(temp_ExptX.^2 , maskS.anatomyTemplate.^2,[0 0]);  %squaring seems to really help sometimes
    
    %Shift mask from other experiment to align to this one
    maskS.bwCell{1} = circshift(maskS.bwCell{1},[round(-mbest) round(-nbest) 0]);
    %maskS.bwCell{1} = circshift_continous2(maskS.bwCell{1},(-nbest),(-mbest));
    
    figure(40)
    hold on
    contour(maskS.bwCell{1},[.5 .5],'b')
        
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in addCell.
function addCell_Callback(hObject, eventdata, handles)
% hObject    handle to addCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS ACQinfo

[xmicperpix ymicperpix] = getImResolution;

if ~isfield(maskS,'bwCell')
    maskS.bwCell{1} = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
end

if get(handles.indicator,'value') == 1  %Manual Polygon
    %imagesc(double(maskS.im{1}),'Buttondownfcn',@OGBROICallback), colormap gray
    
    breakflag = 0;
    while breakflag == 0
        figure(40)
        imagesc(double(maskS.im{1})), colormap gray
        hold on, contour(maskS.bwCell{1},[.5 .5],'r')
        bwdum = roipoly;
        maskS.bwCell{1} = sign(maskS.bwCell{1}+bwdum);
    end
    
elseif get(handles.indicator,'value') == 2  %Manual circle
    
    figure(40)
    imagesc(double(maskS.im{1}),'Buttondownfcn',@circleROICallback), colormap gray
    
else    
    %Bandpass filter
%     im = maskS.im{1};
%     
%     LPsig = 1.5*(1.6/res);
%     HPsig = 2*(1.6/res);
%     
%     LP = fspecial('gaussian',size(im),1); LP = LP/sum(LP(:));
%     HP = fspecial('gaussian',size(im),2); HP = HP/sum(HP(:));
%     BP = LP-HP;
%     imH = ifft2(fft2(BP).*fft2(im));
%     imH = fftshift(fftshift(imH,1),2);
%     maskS.BP = imH;
%     
%     maskS.BP = LocalZ(maskS.BP,20);
%     
    figure(40)
    imagesc(double(maskS.im{1}),'Buttondownfcn',@OGBROICallback), colormap gray
    
end

hold on
contour(double(maskS.bwCell{1}),[.5 .5],'r')
hold off
axis image

% --- Executes on button press in removeCell.
function removeCell_Callback(hObject, eventdata, handles)
% hObject    handle to removeCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

figure(40), 

imagesc(double(maskS.im{1}),'Buttondownfcn',@removeCellCallback), colormap gray
hold on
contour(double(maskS.bwCell{1}),[.5 .5],'r','Buttondownfcn',@removeCellCallback)
hold off
axis image


% --- Executes on button press in pickGlia.
function pickGlia_Callback(hObject, eventdata, handles)
% hObject    handle to pickGlia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

maskS.bwCell{2} = zeros(size(maskS.im{2}));

figure(41), 
imagesc(maskS.im{2}), colormap gray

SelectGlia



% --- Executes on button press in revCorrGUI.
function revCorrGUI_Callback(hObject, eventdata, handles)
% hObject    handle to revCorrGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pRev


function SelectWin_Callback(hObject, eventdata, handles)
% hObject    handle to SelectWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SelectWin as text
%        str2double(get(hObject,'String')) returns contents of SelectWin as a double


% --- Executes during object creation, after setting all properties.
function SelectWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function blockSize_Callback(hObject, eventdata, handles)
% hObject    handle to blockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blockSize as text
%        str2double(get(hObject,'String')) returns contents of blockSize as a double


% --- Executes during object creation, after setting all properties.
function blockSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);

global G_handles G_RChandles ACQinfo Analyzer maskS f0m cellS bw TCWin symbolInfo bcond funcmap

clear global G_handles G_RChandles ACQinfo Analyzer maskS f0m cellS bw TCWin symbolInfo bcond funcmap



% --- Executes on selection change in indicator.
function indicator_Callback(hObject, eventdata, handles)
% hObject    handle to indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns indicator contents as cell array
%        contents{get(hObject,'Value')} returns selected item from indicator


% --- Executes during object creation, after setting all properties.
function indicator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadMaskTemplate.
function loadMaskTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to loadMaskTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS ACQinfo

maskS = struct;

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
    if isfield(S,'axismap')
        axismap = S.axismap;    %f0m is a cell array with images from each condition
        
        
        %Rescale and threshold before applying as a template
        mag = abs(axismap);
        mag = mag.^2;
        
        ma = prctile(mag(:),99.5);
        mag(find(mag>ma)) = ma;
        mi = prctile(mag(:),.5);
        mag(find(mag<mi)) = mi;
        
        mag = mag-min(mag(:));
        mag = mag/max(mag(:));
        
        
        mag = mag-mean(mag(:));
        mag = mag/std(mag(:));
        maskS.imZ{1} = mag;
        maskS.im{1} = mag;
    elseif isfield(S.maskS,'im')
        
        maskS.imZ = S.maskS.im;
        maskS.im = S.maskS.im;
        
    end
    
    
end

if ~isfield(maskS,'bw')
    maskS.bw = ones(size(maskS.imZ{1}));
end

[xmicperpix ymicperpix] = getImResolution;

maskS.bwCell{1} = zeros(size(maskS.imZ{1}));

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
figure(40), 
imagesc(xdom,ydom,maskS.im{1}), colormap gray
hold on
contour(xdom,ydom,maskS.bwCell{1},[.5 .5],'r')
axis image
hold off



% --- Executes on selection change in templateType.
function templateType_Callback(hObject, eventdata, handles)
% hObject    handle to templateType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns templateType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from templateType


% --- Executes during object creation, after setting all properties.
function templateType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to templateType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskTemplateDelim_Callback(hObject, eventdata, handles)
% hObject    handle to maskTemplateDelim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskTemplateDelim as text
%        str2double(get(hObject,'String')) returns contents of maskTemplateDelim as a double


% --- Executes during object creation, after setting all properties.
function maskTemplateDelim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskTemplateDelim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sbxRegistration.
function sbxRegistration_Callback(hObject, eventdata, handles)
% hObject    handle to sbxRegistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sbxaligntool
%choice = questdlg('Would you register this entire experiment and save it to disk?','No thank you','No thank you');

% --- Executes on button press in sbxSegmentation.
function sbxSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to sbxSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sbxsegmenttool

% --- Executes on button press in usePreAligned.
function usePreAligned_Callback(hObject, eventdata, handles)
% hObject    handle to usePreAligned (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usePreAligned



% --- Executes on button press in saveTraces.
function saveTraces_Callback(hObject, eventdata, handles)
% hObject    handle to saveTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cellS Analyzer

AUE = [Analyzer.M.anim '_u' Analyzer.M.unit '_' Analyzer.M.expt];
path = 'c:\';
filename = strcat(path,AUE,'_cellS');
uisave('cellS',filename)

function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on computeF0 and none of its controls.
function computeF0_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to computeF0 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AutoMask.
function AutoMask_Callback(hObject, eventdata, handles)
% hObject    handle to AutoMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


AutoMask


% --- Executes on button press in frameRemoval.
function frameRemoval_Callback(hObject, eventdata, handles)
% hObject    handle to frameRemoval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frameRemoval



function deletionThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to deletionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deletionThreshold as text
%        str2double(get(hObject,'String')) returns contents of deletionThreshold as a double


% --- Executes during object creation, after setting all properties.
function deletionThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deletionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in anatomyTemp.
function anatomyTemp_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

T = 2; %seconds
nT = 15; %N windows

[bestWindow maskS.anatomyTemplate] = findTempWindow(T,nT);

mi = prctile(maskS.anatomyTemplate(:),1);
ma = prctile(maskS.anatomyTemplate(:),99.8);

figure,imagesc(maskS.anatomyTemplate(:,:),[mi ma]), colormap gray
title('Anatomy Template')
axis image


% --- Executes on button press in motionInfo.
function motionInfo_Callback(hObject, eventdata, handles)
% hObject    handle to motionInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


experimentMotionAnalyzer2


% --- Executes on button press in PCremoval.
function PCremoval_Callback(hObject, eventdata, handles)
% hObject    handle to PCremoval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PCremoval



function nPC_Callback(hObject, eventdata, handles)
% hObject    handle to nPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nPC as text
%        str2double(get(hObject,'String')) returns contents of nPC as a double


% --- Executes during object creation, after setting all properties.
function nPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DSflag.
function DSflag_Callback(hObject, eventdata, handles)
% hObject    handle to DSflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DSflag


% --- Executes on button press in OpticFlowCorrection.
function OpticFlowCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to OpticFlowCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OpticFlowCorrection


% --- Executes on button press in cleanMask.
function cleanMask_Callback(hObject, eventdata, handles)
% hObject    handle to cleanMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

%Clean it. Otherwise you get wierd little pixels in there.
% se = strel('disk',1);
% maskS.bwCell{1} = imopen(maskS.bwCell{1},se);

imlabel = bwlabel(maskS.bwCell{1},4);
cellID = unique(imlabel);
cellID = cellID(2:end);

mi = 10;
ma = 300;
nremoval = 0;
for i = 1:length(cellID)
   
    idx = find(imlabel == cellID(i));
    N = length(idx);
    if N > ma | N < mi
        maskS.bwCell{1}(idx) = 0;
        
        nremoval = nremoval+1
        'Removing mask'
    end
    
end

figure
contour(double(maskS.bwCell{1}),[.5 .5],'r')
hold off
axis image
axis ij
