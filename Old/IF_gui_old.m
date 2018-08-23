function varargout = IF_gui(varargin)
% IF_GUI MATLAB code for IF_gui.fig
%      IF_GUI, by itself, creates a new IF_GUI or raises the existing
%      singleton*.
%
%      H = IF_GUI returns the handle to a new IF_GUI or the handle to
%      the existing singleton*.
%
%      IF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IF_GUI.M with the given input arguments.
%
%      IF_GUI('Property','Value',...) creates a new IF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IF_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IF_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IF_gui

% Last Modified by GUIDE v2.5 20-Aug-2018 10:37:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IF_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @IF_gui_OutputFcn, ...
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


% --- Executes just before IF_gui is made visible.
function IF_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IF_gui (see VARARGIN)

% Choose default command line output for IF_gui
handles.output = hObject;



iniFile = fullfile(cd, 'IFgui_ini.mat');
if exist(iniFile, 'file')==2 %If it exists as a file
    initialValues = load(iniFile);
    if isfield(initialValues,'lastUsedInFolder')
        handles.file.path = initialValues.lastUsedInFolder; %From magicgui.m//matlabcentral.com
    else handles.file.path = cd;
    end
    if isfield(initialValues,'lastUsedOutFolder')
        handles.file.outpath = initialValues.lastUsedOutFolder;
    else handles.file.outpath = cd;
    end
else
    handles.file.path    = cd;
    handles.file.outpath = cd;
    lastUsedInFolder     = cd;
    lastUsedOutFolder    = cd;
    save('IFgui_ini.mat', 'lastUsedInFolder', 'lastUsedOutFolder');
end

handles.file.slashtype    = '/';
handles.file.codepath     = cd;
pfinder                   = strfind(handles.file.codepath, handles.file.slashtype);
handles.file.codeparent   = handles.file.codepath(1:pfinder(end));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IF_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IF_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in rdbt_boxPlot.
function rdbt_boxPlot_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_boxPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_boxPlot


% --- Executes on button press in rdbt_scatPlot.
function rdbt_scatPlot_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_scatPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_scatPlot


% --- Executes on selection change in pop_xScat.
function pop_xScat_Callback(hObject, eventdata, handles)
% hObject    handle to pop_xScat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_xScat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_xScat


% --- Executes during object creation, after setting all properties.
function pop_xScat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_xScat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_yScat.
function pop_yScat_Callback(hObject, eventdata, handles)
% hObject    handle to pop_yScat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_yScat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_yScat


% --- Executes during object creation, after setting all properties.
function pop_yScat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_yScat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_xlimLow_Callback(hObject, eventdata, handles)
% hObject    handle to etx_xlimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_xlimLow as text
%        str2double(get(hObject,'String')) returns contents of etx_xlimLow as a double


% --- Executes during object creation, after setting all properties.
function etx_xlimLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_xlimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_xlimHi_Callback(hObject, eventdata, handles)
% hObject    handle to etx_xlimHi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_xlimHi as text
%        str2double(get(hObject,'String')) returns contents of etx_xlimHi as a double


% --- Executes during object creation, after setting all properties.
function etx_xlimHi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_xlimHi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_ylimLow_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ylimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ylimLow as text
%        str2double(get(hObject,'String')) returns contents of etx_ylimLow as a double


% --- Executes during object creation, after setting all properties.
function etx_ylimLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ylimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_ylimHi_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ylimHi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ylimHi as text
%        str2double(get(hObject,'String')) returns contents of etx_ylimHi as a double


% --- Executes during object creation, after setting all properties.
function etx_ylimHi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ylimHi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch1calc.
function pop_ch1calc_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch1calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch1calc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch1calc


% --- Executes during object creation, after setting all properties.
function pop_ch1calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch1calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch2calc.
function pop_ch2calc_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch2calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch2calc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch2calc


% --- Executes during object creation, after setting all properties.
function pop_ch2calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch2calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch3calc.
function pop_ch3calc_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch3calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch3calc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch3calc


% --- Executes during object creation, after setting all properties.
function pop_ch3calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch3calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch4calc.
function pop_ch4calc_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch4calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch4calc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch4calc


% --- Executes during object creation, after setting all properties.
function pop_ch4calc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch4calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_nucmaskPrefix_Callback(hObject, eventdata, handles)
% hObject    handle to etx_nucmaskPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_nucmaskPrefix as text
%        str2double(get(hObject,'String')) returns contents of etx_nucmaskPrefix as a double


% --- Executes during object creation, after setting all properties.
function etx_nucmaskPrefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_nucmaskPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_cytmaskPrefix_Callback(hObject, eventdata, handles)
% hObject    handle to etx_cytmaskPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_cytmaskPrefix as text
%        str2double(get(hObject,'String')) returns contents of etx_cytmaskPrefix as a double


% --- Executes during object creation, after setting all properties.
function etx_cytmaskPrefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_cytmaskPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch1.
function pop_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch1


% --- Executes during object creation, after setting all properties.
function pop_ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch2.
function pop_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch2


% --- Executes during object creation, after setting all properties.
function pop_ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch3.
function pop_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch3


% --- Executes during object creation, after setting all properties.
function pop_ch3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ch4.
function pop_ch4_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch4


% --- Executes during object creation, after setting all properties.
function pop_ch4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%**************************************************************************

function etx_ch1lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch1lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ch1lbl as text
%        str2double(get(hObject,'String')) returns contents of etx_ch1lbl as a double
templist{1} = get(etx_ch1lbl, 'String')

% --- Executes during object creation, after setting all properties.
function etx_ch1lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ch1lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_ch2lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch2lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ch2lbl as text
%        str2double(get(hObject,'String')) returns contents of etx_ch2lbl as a double


% --- Executes during object creation, after setting all properties.
function etx_ch2lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ch2lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_ch3lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch3lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ch3lbl as text
%        str2double(get(hObject,'String')) returns contents of etx_ch3lbl as a double


% --- Executes during object creation, after setting all properties.
function etx_ch3lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ch3lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_ch4lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch4lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_ch4lbl as text
%        str2double(get(hObject,'String')) returns contents of etx_ch4lbl as a double


% --- Executes during object creation, after setting all properties.
function etx_ch4lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_ch4lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in "Load Source Image".
function pushbutton1_Callback(hObject, eventdata, handles)
handles.file.path = uigetdir(handles.file.path, 'Select the experiment''s Main Folder');
if handles.file.path ~= 0      % Assign the value if they didn't click cancel.
    set(handles.list_inputfiles, 'String', {})
    handles.file.path(handles.file.path=='\')='/';
    % Save the image folder in our ini file.
    lastUsedInFolder = handles.file.path;
    save('IFgui_ini.mat', 'lastUsedInFolder', '-append');
    
    ctr = 1;
    treatmentfolder = dir(handles.file.path);
    for aa = 3:length(treatmentfolder)
        tt_finder = strfind(treatmentfolder(aa).name(1:3), 'an_');
        if tt_finder == 1
            handles.file.treatmentfold{ctr} = char(treatmentfolder(aa).name);
            ctr = ctr+1;
        end
    end
end
set(handles.list_inputfiles,'Value',1); %Add this line to set the active item in the list box to the first item before repopulating the listbox. (Throws error otherwise) 
set(handles.list_inputfiles, 'String', handles.file.treatmentfold)
set(handles.etx_sourceFold, 'String', handles.file.path)
guidata(hObject, handles);

% --- Executes on Load button press for Output.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.file.outpath = uigetdir(handles.file.outpath, 'Select the experiment''s Main Folder');
if handles.file.outpath ~= 0
    % Assign the value if they didn't click cancel.
    handles.file.outpath(handles.file.outpath=='\')='/';
    % Save the image folder in our ini file.
    lastUsedOutFolder = handles.file.outpath;
    save('IFgui_ini.mat', 'lastUsedOutFolder', '-append');
end
set(handles.etx_outFold, 'String', lastUsedOutFolder);
guidata(hObject, handles);



function etx_sourceFold_Callback(hObject, eventdata, handles)
% hObject    handle to etx_sourceFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_sourceFold as text
%        str2double(get(hObject,'String')) returns contents of etx_sourceFold as a double


% --- Executes during object creation, after setting all properties.
function etx_sourceFold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_sourceFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etx_outFold_Callback(hObject, eventdata, handles)
% hObject    handle to etx_outFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_outFold as text
%        str2double(get(hObject,'String')) returns contents of etx_outFold as a double


% --- Executes during object creation, after setting all properties.
function etx_outFold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_outFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_xCorr.
function chk_xCorr_Callback(hObject, eventdata, handles)
% hObject    handle to chk_xCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_xCorr


% --- Executes on button press in rdbt_calcCh1.
function rdbt_calcCh1_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_calcCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_calcCh1


% --- Executes on button press in rdbt_calcCh2.
function rdbt_calcCh2_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_calcCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_calcCh2


% --- Executes on button press in rdbt_calcCh3.
function rdbt_calcCh3_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_calcCh3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_calcCh3


% --- Executes on button press in rdbt_calcCh4.
function rdbt_calcCh4_Callback(hObject, eventdata, handles)
% hObject    handle to rdbt_calcCh4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdbt_calcCh4


% --- Executes on button press in chk_singlescat.
function chk_singlescat_Callback(hObject, eventdata, handles)
% hObject    handle to chk_singlescat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_singlescat


% --- Executes on button press in chk_consolidscat.
function chk_consolidscat_Callback(hObject, eventdata, handles)
% hObject    handle to chk_consolidscat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_consolidscat


% --- Executes on button press in pbt_RUN.
function pbt_RUN_Callback(hObject, eventdata, handles)
% hObject    handle to pbt_RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in chk_treatmentnums.
function chk_treatmentnums_Callback(hObject, eventdata, handles)
% hObject    handle to chk_treatmentnums (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_treatmentnums



function etx_treatmentnums_Callback(hObject, eventdata, handles)
% hObject    handle to etx_treatmentnums (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etx_treatmentnums as text
%        str2double(get(hObject,'String')) returns contents of etx_treatmentnums as a double


% --- Executes during object creation, after setting all properties.
function etx_treatmentnums_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etx_treatmentnums (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_inputfiles.
function list_inputfiles_Callback(hObject, eventdata, handles)
% hObject    handle to list_inputfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_inputfiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_inputfiles


% --- Executes during object creation, after setting all properties.
function list_inputfiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_inputfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
