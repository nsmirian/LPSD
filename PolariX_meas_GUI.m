function varargout = PolariX_meas_GUI(varargin)
% POLARIX_MEAS_GUI MATLAB code for PolariX_meas_GUI.fig
%      POLARIX_MEAS_GUI, by itself, creates a new POLARIX_MEAS_GUI or raises the existing
%      singleton*.
%
%      H = POLARIX_MEAS_GUI returns the handle to a new POLARIX_MEAS_GUI or the handle to
%      the existing singleton*.
%
%      POLARIX_MEAS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POLARIX_MEAS_GUI.M with the given input arguments.
%
%      POLARIX_MEAS_GUI('Property','Value',...) creates a new POLARIX_MEAS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PolariX_meas_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PolariX_meas_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PolariX_meas_GUI

% Last Modified by GUIDE v2.5 25-Aug-2023 13:34:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PolariX_meas_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PolariX_meas_GUI_OutputFcn, ...
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


% --- Executes just before PolariX_meas_GUI is made visible.
function PolariX_meas_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PolariX_meas_GUI (see VARARGIN)

% Choose default command line output for PolariX_meas_GUI
handles.output = hObject;

Update_TimeCalInfo(handles)
Update_EnergyCalInfo(handles)

% set the default settings at the opening of the GUI
% handles = initialize_default_settings(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PolariX_meas_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% % Set the default settings at the opening of the GUI
% function handle_new = initialize_default_settings(handle_old)
% % default choices at opening
% handle_new = handle_old;
% 
% % for calibration settings
% 
% chosen = get(handle_old.bgrMethod,'Value')-1;
% handle_new.CalibrationBgrMethod = chosen;  % default background method

% chosen = str2double(get(handle_old.edit_CalibrationNumBGR,'String'));
% handle_new.CalibrationBgrNum = chosen;  % default background number
% 
% chosen = str2double(get(handle_old.edit_CalibrationNumIMG,'String'));
% handle_new.CalibrationImgNum = chosen;  % default image number
% 
% chosen = str2double(get(handle_old.edit_CalibrationScanStep,'String'));
% handle_new.CalibrationScanStep = chosen;  % default scan steps
% 
% chosen = get(handle_old.checkbox_CalibrationWaitForTimeReso,'Value');
% handle_new.CalibrationWaitForTimeReso = chosen;  % default setting for waiting before proceed with time resolution measurement
% 
% % for measurement settings
% chosen = get(handle_old.popupmenu_MeasurementBgrMethod,'Value')-1;
% handle_new.MeasurementBgrMethod = chosen;  % default background method
% 
% chosen = str2double(get(handle_old.edit_MeasurementNumBGR,'String'));
% handle_new.MeasurementBgrNum = chosen;  % default background number
% 
% chosen = str2double(get(handle_old.edit_MeasurementNumIMG,'String'));
% handle_new.MeasurementImgNum = chosen;  % default image number
% 
% % for measurement settings
% contents = cellstr(get(handle_old.popupmenu_gainMode,'String'));
% gainMode = contents{get(handle_old.popupmenu_gainMode,'Value')};
% handle_new.gainMode = gainMode;  % default adjust gain mode


% --- Outputs from this function are returned to the command line.
function varargout = PolariX_meas_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in doTimeCali.
function doTimeCali_Callback(hObject, eventdata, handles)
% hObject    handle to doTimeCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.phaseLeft,'String')) && ~isempty(get(handles.phaseRight,'String'))

    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING TIME CALIBRATION - sit back and relax. :-)');
    
	phase_left = str2double(get(handles.phaseLeft,'String'));
	phase_right = str2double(get(handles.phaseRight,'String'));
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
    numStep = str2double(get(handles.numStep,'String'));
    
    comment = get(handles.comment,'String');

    try
        get_timecalibration_OTR9FL2XTDS(phase_left, phase_right, bgrMethod, numBgr, numSig, numStep, msgLoc, comment)
    end
    
    Update_TimeCalInfo(handles)
    display_message(msgLoc, 'FINISHED TIME CALIBRATION.');
end


% --- Executes on button press in doTimeCaliFast.
function doTimeCaliFast_Callback(hObject, eventdata, handles)
% hObject    handle to doTimeCaliFast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.phaseLeft,'String')) && ~isempty(get(handles.phaseRight,'String'))

    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING FAST TIME CALIBRATION  - sit back and relax. :-)');
    
	phase_left = str2double(get(handles.phaseLeft,'String'));
	phase_right = str2double(get(handles.phaseRight,'String'));
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
    numStep = str2double(get(handles.numStep,'String'));
    
    comment = get(handles.comment,'String');
    
    try
        get_timecalibration_fast_OTR9FL2XTDS(phase_left, phase_right, bgrMethod, numBgr, numSig, numStep, msgLoc, comment)
    end
    
    Update_TimeCalInfo(handles)
    display_message(msgLoc, 'FINISHED FAST TIME CALIBRATION.');
end

% --- Executes on button press in doTimeCaliFastDouble.
function doTimeCaliFastDouble_Callback(hObject, eventdata, handles)
% hObject    handle to doTimeCaliFastDouble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.phaseLeft,'String')) && ~isempty(get(handles.phaseRight,'String'))

    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING FAST TIME CALIBRATION in two directions - sit back and relax. :-)');
    
	phase_left = str2double(get(handles.phaseLeft,'String'));
	phase_right = str2double(get(handles.phaseRight,'String'));
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
    numStep = str2double(get(handles.numStep,'String'));
    
    comment = get(handles.comment,'String');
    
    try
        get_timecalibration_fast_double_OTR9FL2XTDS(phase_left, phase_right, bgrMethod, numBgr, numSig, numStep, msgLoc, comment)
    end
    
    Update_TimeCalInfo(handles)
    
    display_message(msgLoc, 'FINISHED FAST TIME CALIBRATION in two directions.');
end

% --- Executes on button press in doEnergyCali.
function doEnergyCali_Callback(hObject, eventdata, handles)
% hObject    handle to doEnergyCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.currentTop,'String')) && ~isempty(get(handles.currentBottom,'String'))

    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING ENERGY CALIBRATION - sit back and relax. :-)');
    
	currentTop = str2double(get(handles.currentTop,'String'));
	currentBottom = str2double(get(handles.currentBottom,'String'));
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
    numStep = str2double(get(handles.numStep,'String'));
    
    comment = get(handles.comment,'String');
    
    try
        get_energycalibration_OTR9FL2XTDS(currentTop, currentBottom, bgrMethod, numBgr, numSig, numStep, msgLoc, comment)
    end
    
    Update_EnergyCalInfo(handles)
    
    display_message(msgLoc, 'FINISHED ENERGY CALIBRATION.');
end
%%%%%%
% we need to add another collibration here 
function energy_colli_with_cycling_Callback(hObject, eventdata, handles)
% hObject    handle to doEnergyCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(get(handles.currentTop,'String')) && ~isempty(get(handles.currentBottom,'String'))

    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING ENERGY CALIBRATION - sit back and relax. :-)');
    
	currentTop = str2double(get(handles.currentTop,'String'));
	currentBottom = str2double(get(handles.currentBottom,'String'));
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
    numStep = str2double(get(handles.numStep,'String'));
     
    comment = get(handles.comment,'String');
    
    try
        get_energycalibration_cycling(currentTop, currentBottom, bgrMethod, numBgr, numSig, numStep, msgLoc, comment)
    end
    
    Update_EnergyCalInfo(handles)
    
    display_message(msgLoc, 'FINISHED ENERGY CALIBRATION.');
end
%%%%%
% --- Executes on button press in doMeasurement.
function doMeasurement_Callback(hObject, eventdata, handles)
% hObject    handle to doMeasurement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    msgLoc = handles.txtOut;
    display_message(msgLoc, 'STARTING DATA TAKING - sit back and relax. :-)');
    
    bgrMethod = get(handles.bgrMethod,'Value')-1;
    numBgr = str2double(get(handles.numBgr,'String'));
    numSig = str2double(get(handles.numSig,'String'));
%     numStep = str2double(get(handles.numStep,'String'));
    
    comment = get(handles.comment,'String');
    get_calibratedimage_OTR9FL2XTDS(bgrMethod, numBgr, numSig, msgLoc, comment)  
    display_message(msgLoc, 'FINISHED DATA TAKING');


% --- Executes on button press in loadTimeCali.
function loadTimeCali_Callback(hObject, eventdata, handles)
% hObject    handle to loadTimeCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

strPattern = ['timecalibration*.mat'];
dialogTitle = 'time calibration file';
folderpath = ['time_calibration/'];

filename_timecal = selectFiles_for_popupmenu(hObject, strPattern,dialogTitle,folderpath);
load([folderpath, filename_timecal], 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'timestamp')
save('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FLFXTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'timestamp')
Update_TimeCalInfo(handles)

% --- Executes on button press in loadEnergyCali.
function loadEnergyCali_Callback(hObject, eventdata, handles)
% hObject    handle to loadEnergyCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
strPattern = ['energycalibration*.mat'];
dialogTitle = 'energy calibration file';
folderpath = ['energy_calibration/'];

filename_energycal = selectFiles_for_popupmenu(hObject, strPattern,dialogTitle,folderpath);
load([folderpath, filename_energycal], 'ergcal', 'ergcal_err', 'timestamp')
save('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/erg_calib_9FLFXTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'timestamp')
Update_TimeCalInfo(handles)

function filename = selectFiles_for_popupmenu(hPopupmenu, strPattern,dialogTitle,folderpath)

% filename_old = get(hPopupmenu,'String');

[filename,pathname] = uigetfile(strPattern, dialogTitle, folderpath);


function phaseRight_Callback(hObject, eventdata, handles)
% hObject    handle to phaseRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseRight as text
%        str2double(get(hObject,'String')) returns contents of phaseRight as a double


% --- Executes during object creation, after setting all properties.
function phaseRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phaseLeft_Callback(hObject, eventdata, handles)
% hObject    handle to phaseLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseLeft as text
%        str2double(get(hObject,'String')) returns contents of phaseLeft as a double


% --- Executes during object creation, after setting all properties.
function phaseLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in refreshTimeCali.
function refreshTimeCali_Callback(hObject, eventdata, handles)
% hObject    handle to refreshTimeCali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% display_message(handles.txtOut , 'test')

Update_TimeCalInfo(handles)


% --- Executes on button press in readLeftPhase.
function readLeftPhase_Callback(hObject, eventdata, handles)
% hObject    handle to readLeftPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% read left phase
try
    PolariX_phase = 'FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
    current_phase = hlcr(PolariX_phase);
    set(handles.phaseLeft,'String',num2str(current_phase, '%.1f'));
end


% --- Executes on button press in readRightPhase.
function readRightPhase_Callback(hObject, eventdata, handles)
% hObject    handle to readRightPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% read left phase
try
    PolariX_phase = 'FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
    current_phase = hlcr(PolariX_phase);
    set(handles.phaseRight,'String',num2str(current_phase, '%.1f'));
end




function currentBottom_Callback(hObject, eventdata, handles)
% hObject    handle to currentBottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentBottom as text
%        str2double(get(hObject,'String')) returns contents of currentBottom as a double


% --- Executes during object creation, after setting all properties.
function currentBottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentBottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentTop_Callback(hObject, eventdata, handles)
% hObject    handle to currentTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentTop as text
%        str2double(get(hObject,'String')) returns contents of currentTop as a double


% --- Executes during object creation, after setting all properties.
function currentTop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in readTopCurrent.
function readTopCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to readTopCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    
    %current_current = hlcr('TTF2.MAGNETS/DIPOLE/D9SMATCH/PS');
    current_current = hlcr('FLASH.MAGNETS/MAGNET.ML/D4FL2XTDS/CURRENT.SP');
    set(handles.currentTop,'String',num2str(current_current, '%.1f'));
end

% --- Executes on button press in readBottomCurrent.
function readBottomCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to readBottomCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    
    %current_current = hlcr('TTF2.MAGNETS/DIPOLE/D9SMATCH/PS');
    current_current = hlcr('FLASH.MAGNETS/MAGNET.ML/D4FL2XTDS/CURRENT.SP');
    set(handles.currentBottom,'String',num2str(current_current, '%.1f'));
end




function numStep_Callback(hObject, eventdata, handles)                                                                       % Not clear to me (N)
% hObject    handle to numStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numStep as text
%        str2double(get(hObject,'String')) returns contents of numStep as a double


% --- Executes during object creation, after setting all properties.
function numStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numSig_Callback(hObject, eventdata, handles)                                                                       % Not clear to me (N)
% hObject    handle to numSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numSig as text
%        str2double(get(hObject,'String')) returns contents of numSig as a double


% --- Executes during object creation, after setting all properties.
function numSig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numBgr_Callback(hObject, eventdata, handles)
% hObject    handle to numBgr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numBgr as text
%        str2double(get(hObject,'String')) returns contents of numBgr as a double


% --- Executes during object creation, after setting all properties.
function numBgr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numBgr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bgrMethod.
function bgrMethod_Callback(hObject, eventdata, handles)
% hObject    handle to bgrMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bgrMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bgrMethod


% --- Executes during object creation, after setting all properties.
function bgrMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bgrMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Update_TimeCalInfo(handles)
try
    load('/home/ttflinac/user/mflorian/PolariX/FL2/time_calibration/time_calib_9FL2XTDS.mat', 'amplitude_XTDS', 'timecal_fspixel', 'timecal_fspixel_err', 'timecal_streak', 'timestamp')

    set(handles.timeCaliTimeStamp,'String',timestamp);
    set(handles.timeCaliValue,'String',[num2str(timecal_fspixel, '%.3f'), '+-', num2str(timecal_fspixel_err, '%.3f') ' fs/pixel, Streak: ', num2str(timecal_streak, '%.1f')])
catch
    set(handles.timeCaliTimeStamp,'String','no old time cali file available');
end

function Update_EnergyCalInfo(handles)
try
    load('/home/ttflinac/user/mflorian/PolariX/FL2/energy_calibration/erg_calib_9FL2XTDS.mat', 'ergcal', 'ergcal_err', 'timestamp')
catch
    set(handles.timeCaliTimeStamp,'String','no old energy cali file available');
end

set(handles.energyCaliTimeStamp,'String',timestamp);
set(handles.energyCaliValue,'String',[num2str(ergcal*10, '%.3f'), '+-', num2str(ergcal_err*10, '%.3f') ' *1e-4/pixel'])


% --- Executes on selection change in txtOut.
function txtOut_Callback(hObject, eventdata, handles)
% hObject    handle to txtOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns txtOut contents as cell array
%        contents{get(hObject,'Value')} returns selected item from txtOut


% --- Executes during object creation, after setting all properties.
function txtOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function comment_Callback(hObject, eventdata, handles)
% hObject    handle to comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comment as text
%        str2double(get(hObject,'String')) returns contents of comment as a double


% --- Executes during object creation, after setting all properties.
function comment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
