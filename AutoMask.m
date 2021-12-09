function varargout = AutoMask(varargin)
% AUTOMASK MATLAB code for AutoMask.fig
%      AUTOMASK, by itself, creates a new AUTOMASK or raises the existing
%      singleton*.
%
%      H = AUTOMASK returns the handle to a new AUTOMASK or the handle to
%      the existing singleton*.
%
%      AUTOMASK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOMASK.M with the given input arguments.
%
%      AUTOMASK('Property','Value',...) creates a new AUTOMASK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoMask_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoMask_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoMask

% Last Modified by GUIDE v2.5 23-Nov-2017 14:04:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoMask_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoMask_OutputFcn, ...
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


% --- Executes just before AutoMask is made visible.
function AutoMask_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoMask (see VARARGIN)

% Choose default command line output for AutoMask
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoMask wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AutoMask_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function maxMinSize_Callback(hObject, eventdata, handles)
% hObject    handle to maxMinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxMinSize as text
%        str2double(get(hObject,'String')) returns contents of maxMinSize as a double


% --- Executes during object creation, after setting all properties.
function maxMinSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxMinSize (see GCBO)
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


% --- Executes on button press in makeMask.
function makeMask_Callback(hObject, eventdata, handles)
% hObject    handle to makeMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS ACQinfo

%'f0forMask' uses the entire experiment to compute the template for the mask.
%N.B. Process.raw must have already been run to compute the 'f0m' images...
%and w/o baseline subtraction

[xmicperpix ymicperpix] = getImResolution;

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;

imZ = Znormalize(str2num(get(handles.maskSize,'string'))); %Take local Zscore

imZ = imZ - ones(size(imZ,1),1)*mean(imZ);

maskMorph = 3;
%maskMorph = str2num(get(handles.maskMorph,'string'));

maskS.bwCell{1} = ZThresh(imZ,str2num(get(handles.maskThresh,'string')),maskMorph);

%minimum cell size (should be done after applying ROI)
eval(['mima = ' get(handles.maxMinSize,'string') ' '])
maskS.bwCell{1} = cellMinMaxSize(maskS.bwCell{1},mima);
maskS.bwCell{1} = cellMorph(maskS.bwCell{1},maskMorph); 

%minimum cell size (should be done after applying ROI)
eval(['mima = ' get(handles.maxMinSize,'string') ' '])
maskS.bwCell{1} = cellMinMaxSize(maskS.bwCell{1},mima);

hold on
contour(maskS.bwCell{1},[.5 .5],'r')

figure,imagesc(xdom,ydom,maskS.bwCell{1}), colormap gray

button = questdlg('Reset region of interest?','ROI','Yes');
if strcmp(button,'Yes')
    dum = roipoly;
    maskS.bwCell{1} = maskS.bwCell{1}.*dum;
end

mi = prctile(maskS.im{1}(:),0);
ma = prctile(maskS.im{1}(:),100);
dum = maskS.im{1};
dum = (dum-mi)/(ma-mi);
dum(find(dum>1)) = 1;
dum(find(dum<0)) = 0;

figure(80), 
imagesc(xdom,ydom,dum), colormap gray
hold on
contour(xdom,ydom,maskS.bwCell{1},[.5 .5],'r')
axis image
hold off


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
