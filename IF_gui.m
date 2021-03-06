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

% Last Modified by GUIDE v2.5 21-Aug-2018 17:49:47

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

set(handles.chk_nucMaskPrefix, 'Value', 1);
set(handles.etx_ch1lbl, 'Enable', 'off');
set(handles.etx_ch2lbl, 'Enable', 'off');
set(handles.etx_ch3lbl, 'Enable', 'off');
set(handles.etx_ch4lbl, 'Enable', 'off');
set(handles.rdbt_noNorm, 'Value', 1);
handles.calcs.normFlag    = 0;

%Add tools
addpath([handles.file.codeparent handles.file.slashtype 'DhruvTools']);
addpath([handles.file.codeparent handles.file.slashtype 'ExchangeTools' handles.file.slashtype 'UnivarScatter_v1.1']);
addpath([handles.file.codeparent handles.file.slashtype 'ExchangeTools' handles.file.slashtype 'DrosteEffect-BrewerMap-a77e675']);
addpath([handles.file.codeparent handles.file.slashtype 'ExchangeTools' handles.file.slashtype 'raacampbell-notBoxPlot-2fbf98c' handles.file.slashtype 'code']);


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




%**************************************************************************

function etx_ch1lbl_Callback(hObject, eventdata, handles)
% Functionality to update the Axes Property popup dialog
%list with all current active channels

if get(handles.pop_ch1, 'Value')>1
    templist{1} = get(handles.etx_ch1lbl, 'String');
else
    templist{1} = [];
end

if get(handles.pop_ch2, 'Value')>1
    templist{2} = get(handles.etx_ch2lbl, 'String');
else
    templist{2} = [];
end

if get(handles.pop_ch3, 'Value')>1
    templist{3} = get(handles.etx_ch3lbl, 'String');
else
    templist{3} = [];
end

if get(handles.pop_ch4, 'Value')>1
    templist{4} = get(handles.etx_ch4lbl, 'String');
else
    templist{4} = [];
end
set(handles.pop_yScat, 'String', templist);
set(handles.pop_xScat, 'String', templist);
set(handles.tx_chlabel1, 'String', templist{1});
guidata(hObject, handles)


function etx_ch2lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch2lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.pop_ch1, 'Value')>1
    templist{1} = get(handles.etx_ch1lbl, 'String');
else
    templist{1} = [];
end

if get(handles.pop_ch2, 'Value')>1
    templist{2} = get(handles.etx_ch2lbl, 'String');
else
    templist{2} = [];
end

if get(handles.pop_ch3, 'Value')>1
    templist{3} = get(handles.etx_ch3lbl, 'String');
else
    templist{3} = [];
end

if get(handles.pop_ch4, 'Value')>1
    templist{4} = get(handles.etx_ch4lbl, 'String');
else
    templist{4} = [];
end
set(handles.pop_yScat, 'String', templist);
set(handles.pop_xScat, 'String', templist);
set(handles.tx_chlabel2, 'String', templist{2});
guidata(hObject, handles)



function etx_ch3lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch3lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.pop_ch1, 'Value')>1
    templist{1} = get(handles.etx_ch1lbl, 'String');
else
    templist{1} = [];
end

if get(handles.pop_ch2, 'Value')>1
    templist{2} = get(handles.etx_ch2lbl, 'String');
else
    templist{2} = [];
end

if get(handles.pop_ch3, 'Value')>1
    templist{3} = get(handles.etx_ch3lbl, 'String');
else
    templist{3} = [];
end

if get(handles.pop_ch4, 'Value')>1
    templist{4} = get(handles.etx_ch4lbl, 'String');
else
    templist{4} = [];
end
set(handles.pop_yScat, 'String', templist);
set(handles.pop_xScat, 'String', templist);
set(handles.tx_chlabel3, 'String', templist{3});
guidata(hObject, handles)


function etx_ch4lbl_Callback(hObject, eventdata, handles)
% hObject    handle to etx_ch4lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.pop_ch1, 'Value')>1
    templist{1} = get(handles.etx_ch1lbl, 'String');
else
    templist{1} = [];
end

if get(handles.pop_ch2, 'Value')>1
    templist{2} = get(handles.etx_ch2lbl, 'String');
else
    templist{2} = [];
end

if get(handles.pop_ch3, 'Value')>1
    templist{3} = get(handles.etx_ch3lbl, 'String');
else
    templist{3} = [];
end

if get(handles.pop_ch4, 'Value')>1
    templist{4} = get(handles.etx_ch4lbl, 'String');
else
    templist{4} = [];
end
set(handles.pop_yScat, 'String', templist);
set(handles.pop_xScat, 'String', templist);
set(handles.tx_chlabel4, 'String', templist{4});
guidata(hObject, handles)


% --- Executes on button press in "Load Source Data".
function pushbutton1_Callback(hObject, eventdata, handles)
handles.file.path = uigetdir(handles.file.path, 'Select the experiment''s Main Folder');

if handles.file.path ~= 0      % Assign the value if they didn't click cancel.
    set(handles.list_inputfiles, 'String', {})
    handles.file.path(handles.file.path=='\')='/';
    
    % Save the image folder in our ini file.
    lastUsedInFolder = handles.file.path;
    save('IFgui_ini.mat', 'lastUsedInFolder', '-append');

    ctr             = 1;
    treatmentfolder = dir(handles.file.path);
    
    for aa = 3:length(treatmentfolder)
        tt_finder = strfind(treatmentfolder(aa).name(1:3), 'an_');
        
        if tt_finder == 1
            handles.file.treatmentfold{ctr} = char(treatmentfolder(aa).name);
            ctr                             = ctr+1;
        end
        
    end
    
    
end


set(handles.list_inputfiles,'Value',1); %Add this line to set the active item in the list box to the first item before repopulating the listbox. (Throws error otherwise)
set(handles.list_inputfiles, 'String', handles.file.treatmentfold)
set(handles.etx_sourceFold, 'String', handles.file.path)
namefinder                  = strfind(handles.file.path, handles.file.slashtype);
handles.file.experimentName = handles.file.path(namefinder(end)+1:end);

%Load metadata:
metaname = dir(fullfile(handles.file.path, '*meta.mat'));
if ~isempty(metaname)
    metafile = ([handles.file.path handles.file.slashtype metaname(1).name]);
    button   = questdlg(['Metadata file found! It''s name is: ' metaname(1).name '  Should I load these settings? '], 'Experiment Metadata Loader', 'Yes', 'No', 'Yes');
    
    if strcmp(button, 'Yes')
        close(gcf)
        load(metafile)
%         copyFields = {'outputs', 'output', 'inputs', 'calcs', 'file'};
%         handles    = LoadGuiState(metafile, handles, copyFields); %no need to run LoadGuiState!! Figured out an easier way.
        
    end
    
end

set(handles.pop_treatList, 'String', cleanNames(handles.file.treatmentfold))
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



% --- Executes on selection change in pop_ch1.
function pop_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch1

if get(handles.pop_ch1, 'Value')~=1
    set(handles.etx_ch1lbl, 'Enable', 'on');
else
    set(handles.etx_ch1lbl, 'Enable', 'off');
end
guidata(hObject, handles)


% --- Executes on selection change in pop_ch2.
function pop_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch2
if get(handles.pop_ch2, 'Value')~=1
    set(handles.etx_ch2lbl, 'Enable', 'on');
else
    set(handles.etx_ch2lbl, 'Enable', 'off');
end
guidata(hObject, handles)



% --- Executes on selection change in pop_ch3.
function pop_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch3
if get(handles.pop_ch3, 'Value')~=1
    set(handles.etx_ch3lbl, 'Enable', 'on');
else
    set(handles.etx_ch3lbl, 'Enable', 'off');
end
guidata(hObject, handles)




% --- Executes on selection change in pop_ch4.
function pop_ch4_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ch4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ch4
if get(handles.pop_ch4, 'Value')~=1
    set(handles.etx_ch4lbl, 'Enable', 'on');
else
    set(handles.etx_ch4lbl, 'Enable', 'off');
end
guidata(hObject, handles)



% --- Executes on button press in pbt_xlimDefault.
function pbt_xlimDefault_Callback(hObject, eventdata, handles)

set(handles.etx_xlimLow, 'String', num2str(0));
set(handles.etx_xlimHi, 'String', num2str(250));
guidata(hObject, handles)

% --- Executes on button press in pbt_ylimDefault.
function pbt_ylimDefault_Callback(hObject, eventdata, handles)

set(handles.etx_ylimLow, 'String', num2str(0));
set(handles.etx_ylimHi, 'String', num2str(250));
guidata(hObject, handles)

%=========================================================================

%*************************** RUN *****************************************

%=========================================================================

% --- Executes on button press in pbt_RUN.
function pbt_RUN_Callback(hObject, eventdata, handles)
% hObject    handle to pbt_RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%File
handles.file.treatmentLabels = cleanNames(handles.file.treatmentfold, '_');

%Inputs
handles.inputs.ch1lab        = get(handles.etx_ch1lbl, 'String');
handles.inputs.ch2lab        = get(handles.etx_ch2lbl, 'String');
handles.inputs.ch3lab        = get(handles.etx_ch3lbl, 'String');
handles.inputs.ch4lab        = get(handles.etx_ch4lbl, 'String');
handles.inputs.nucMaskPrefix = get(handles.etx_nucmaskPrefix, 'String');
handles.inputs.cytMaskPrefix = get(handles.etx_cytmaskPrefix, 'String');
handles.inputs.ChannelLabel  = get(handles.pop_xScat, 'String');
handles.inputs.activeChans   = [get(handles.rdbt_calcCh1, 'Value') ...
    get(handles.rdbt_calcCh2, 'Value') ...
    get(handles.rdbt_calcCh3, 'Value') ...
    get(handles.rdbt_calcCh4, 'Value')];

%Outputs

%ImageFormat
if get(handles.rdbt_png, 'Value')==1
    handles.outputs.imageFormat = 'png';
elseif get(handles.rdbt_svg, 'Value')==1
    handles.outputs.imageFormat = 'svg';
end

%PlotMode
if get(handles.rdbt_plotAll, 'Value')
    handles.outputs.plotMode = 'all';
elseif get(handles.rdbt_plotOnly, 'Value')
    handles.outputs.plotMode = 'subset';
    
    %error handling
    if strcmp(get(handles.etx_treatmentnums, 'String'), '...')
        errordlg('No treatment numbers entered! ')
    end
    
end

handles.outputs.boxAutoY     = get(handles.pbt_boxAutoLim, 'Value');
handles.outputs.boxYlim      = [str2double(get(handles.etx_boxLimLow, 'String')) str2double(get(handles.etx_boxLimHigh, 'String'))];
handles.outputs.boxplot      = get(handles.rdbt_boxPlot, 'Value');
handles.outputs.conscat      = get(handles.chk_consolidscat, 'Value');
handles.outputs.chanscat     = get(handles.chk_singlescat,'Value');
handles.outputs.limconscatID = str2num(get(handles.etx_treatmentnums, 'String')); %#ok<ST2NM>
handles.outputs.scatX        = get(handles.pop_xScat,'Value' );
handles.outputs.scatY        = get(handles.pop_yScat, 'Value');
handles.outputs.scatAutoX    = get(handles.pbt_xlimDefault, 'Value');
handles.outputs.scatAutoY    = get(handles.pbt_ylimDefault, 'Value');
handles.outputs.scatXlim     = [str2double(get(handles.etx_xlimLow, 'String')) str2double(get(handles.etx_xlimHi, 'String'))];
handles.outputs.scatYlim     = [str2double(get(handles.etx_ylimLow, 'String')) str2double(get(handles.etx_ylimHi, 'String'))];
handles.outputs.CorrLine     = get(handles.chk_xCorr, 'Value');
handles.outputs.margDist     = get(handles.chk_margDist, 'Value');

%Calcs
calclist(1:4) = 0;
calclabels = get(handles.pop_ch1calc, 'String');
[handles.calcs.nuc, handles.calcs.cyt, ...
    handles.calcs.wholeCell, ...
    handles.calcs.ncRatio, ...
    handles.calcs.cnRatio]= deal(zeros([1 4]));

if get(handles.rdbt_calcCh1,'Value')
    calclist(1) = get(handles.pop_ch1calc, 'Value');
    handles.calcs.label{1}  = calclabels{calclist(1)};
end

if get(handles.rdbt_calcCh2,'Value')
    calclist(2) = get(handles.pop_ch2calc, 'Value');
    handles.calcs.label{2}  = calclabels{calclist(2)};
end

if get(handles.rdbt_calcCh3,'Value')
    calclist(3) = get(handles.pop_ch3calc, 'Value');
    handles.calcs.label{3}  = calclabels{calclist(3)};
end

if get(handles.rdbt_calcCh4,'Value')
    calclist(4) = get(handles.pop_ch4calc, 'Value');
    handles.calcs.label{4}  = calclabels{calclist(4)};
end

calclist = [calclist];

handles.calcs.nuc      (find(calclist==1)) = 1;  %This can be re-written to be quicker!
handles.calcs.cyt      (find(calclist==2)) = 1;
handles.calcs.wholeCell(find(calclist==3)) = 1;
handles.calcs.ncRatio  (find(calclist==4)) = 1;
handles.calcs.cnRatio  (find(calclist==5)) = 1;

%Edit this variable if you want multiple calculations in the same channel:
handles.calcs.all = [handles.calcs.nuc; handles.calcs.cyt; ...
    handles.calcs.wholeCell; handles.calcs.ncRatio; ...
    handles.calcs.cnRatio];

%Normalization flags:
if get(handles.rdbt_yesNorm, 'Value')
    t_normList            = get(handles.pop_normType, 'String');
   handles.calcs.normType = char(t_normList(get(handles.pop_normType, 'Value'))); 
   handles.calcs.normTo   = get(handles.pop_treatList, 'Value');
   handles.calcs.normFlag = 1;
end

%Save metadata:
handlesLoad = handles;
save([handles.file.path handles.file.slashtype char(cleanNames({handles.file.experimentName}, '_')) 'meta.mat'], 'handlesLoad');


%Call Main Fx
IF_gui_backend(handles.file, handles.inputs, handles.calcs, handles.outputs)


%=========================================================================

%****************** UNUSED CALLBACKS *************************************

%=========================================================================

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)

% --- Executes on button press in chk_treatmentnums.
function chk_treatmentnums_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_boxPlot.
function rdbt_boxPlot_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_scatPlot.
function rdbt_scatPlot_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_xScat.
function pop_xScat_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_yScat.
function pop_yScat_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_ch1calc.
function pop_ch1calc_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_ch2calc.
function pop_ch2calc_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_ch3calc.
function pop_ch3calc_Callback(hObject, eventdata, handles)

% --- Executes on selection change in pop_ch4calc.
function pop_ch4calc_Callback(hObject, eventdata, handles)

% --- Executes on button press in chk_xCorr.
function chk_xCorr_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_calcCh1.
function rdbt_calcCh1_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_calcCh2.
function rdbt_calcCh2_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_calcCh3.
function rdbt_calcCh3_Callback(hObject, eventdata, handles)

% --- Executes on button press in rdbt_calcCh4.
function rdbt_calcCh4_Callback(hObject, eventdata, handles)

% --- Executes on button press in chk_singlescat.
function chk_singlescat_Callback(hObject, eventdata, handles)

% --- Executes on button press in chk_consolidscat.
function chk_consolidscat_Callback(hObject, eventdata, handles)


%-------------------------------------------------------------------------
%-------------- OBJECT CREATION CALLBACKS --------------------------------
%-------------------------------------------------------------------------
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

% --- Textbox Callbacks.
function etx_xlimHi_Callback(hObject, eventdata, handles)

function etx_xlimLow_Callback(hObject, eventdata, handles)

function etx_ylimHi_Callback(hObject, eventdata, handles)

function etx_ylimLow_Callback(hObject, eventdata, handles)

function etx_treatmentnums_Callback(hObject, eventdata, handles)

function etx_nucmaskPrefix_Callback(hObject, eventdata, handles)

function etx_cytmaskPrefix_Callback(hObject, eventdata, handles)
