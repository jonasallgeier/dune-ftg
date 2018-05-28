function varargout = graphical_analysis(varargin)
% GRAPHICAL_ANALYSIS M-file for graphical_analysis.fig
%      GRAPHICAL_ANALYSIS, by itself, creates a new GRAPHICAL_ANALYSIS or raises the existing
%      singleton*.
%
%      H = GRAPHICAL_ANALYSIS returns the handle to a new GRAPHICAL_ANALYSIS or the handle to
%      the existing singleton*.
%
%      GRAPHICAL_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAPHICAL_ANALYSIS.M with the given input arguments.
%
%      GRAPHICAL_ANALYSIS('Property','Value',...) creates a new GRAPHICAL_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simple_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to graphical_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help graphical_analysis

% Copyright 2001-2003 The MathWorks, Inc.

% Last Modified by GUIDE v2.5 18-May-2018 13:51:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @graphical_analysis_OpeningFcn, ...
    'gui_OutputFcn',  @graphical_analysis_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before graphical_analysis is made visible.
function graphical_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to graphical_analysis (see VARARGIN)
% Create the data to plot

% needed for print button
addpath('./export_fig');

[switches,Cancelled] = input_dialog_1;
if Cancelled
    error('User terminated.');
end
[user_input,Cancelled] = input_dialog_2(switches);
if Cancelled
    error('User terminated.');
end

%prompt = {'Number of electrodes','File of electrode coordinates','Number of combinations','Name of electrode configuration file','Results directory/basename'};
%title = 'Input';
%dims = [1 60];
%definput = {'144','electrode.configuration','4650','4_electrodes.txt','./data/longtestrun'};
%answer = inputdlg(prompt,title,dims,definput);


handles.n_conf = str2num(user_input.n_comb);
[~,handles.A_elec,handles.B_elec,handles.M_elec,handles.N_elec]=textread(user_input.file_combs,'%f %f %f %f %f',...
    handles.n_conf,'headerlines',1);

[~,handles.X_elec,handles.Y_elec,handles.Z_elec,~]=textread(user_input.file_coords,'%f %f %f %f %d',...
    str2num(user_input.n_el),'headerlines',1);

handles = read_potentials_transient(user_input,handles);
if switches.doPotentialsBase
    handles = read_potentials_base(user_input,handles);
end
if switches.doPotentialsMoments
   handles = read_potentials_moments(user_input,handles);
end

handles.switches = switches;
% Choose default command line output for graphical_analysis
handles.output = hObject;

set(handles.conf_slider,'Min',1);
set(handles.conf_slider,'Max',handles.n_conf);
set(handles.conf_slider,'Value',1);
set(handles.conf_slider, 'SliderStep', [1/(handles.n_conf-1), 100/(handles.n_conf-1)]);
set(handles.conf_number,'String','1');

set(handles.A_electrode,'string',num2str(handles.A_elec(1),'%d'));
set(handles.M_electrode,'string',num2str(handles.M_elec(1),'%d'));
set(handles.N_electrode,'string',num2str(handles.N_elec(1),'%d'));
set(handles.B_electrode,'string',num2str(handles.B_elec(1),'%d'));

update_figure(handles);

guidata(hObject, handles);
% UIWAIT makes graphical_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = graphical_analysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end




% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end



function A_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_electrode as text
%        str2double(get(hObject,'String')) returns contents of A_electrode as a double
update_figure(handles);
set(handles.conf_number,'String','manual');
end

% --- Executes during object creation, after setting all properties.
function A_electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function M_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to M_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_electrode as text
%        str2double(get(hObject,'String')) returns contents of M_electrode as a double
update_figure(handles);
set(handles.conf_number,'String','manual');
end

% --- Executes during object creation, after setting all properties.
function M_electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function N_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to N_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_electrode as text
%        str2double(get(hObject,'String')) returns contents of N_electrode as a double
update_figure(handles);
set(handles.conf_number,'String','manual');
end

% --- Executes during object creation, after setting all properties.
function N_electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function B_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to B_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_electrode as text
%        str2double(get(hObject,'String')) returns contents of B_electrode as a double
update_figure(handles);
set(handles.conf_number,'String','manual');
end

% --- Executes during object creation, after setting all properties.
function B_electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in apply_changes.
function apply_changes_Callback(hObject, eventdata, handles)
% hObject    handle to apply_changes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_figure(handles);
set(handles.conf_number,'String','manual')
end


% % --- Executes on slider movement.
% function slider1_Callback(hObject, eventdata, handles)
% % hObject    handle to conf_slider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% end

% --- Executes during object creation, after setting all properties.
function conf_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conf_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end


function conf_number_Callback(hObject, eventdata, handles)
% hObject    handle to conf_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conf_number as text
%        str2double(get(hObject,'String')) returns contents of conf_number as a double
val = str2num(get(hObject,'String'));
set(handles.A_electrode,'String',num2str(handles.A_elec(val),'%d'));
set(handles.M_electrode,'String',num2str(handles.M_elec(val),'%d'));
set(handles.N_electrode,'String',num2str(handles.N_elec(val),'%d'));
set(handles.B_electrode,'String',num2str(handles.B_elec(val),'%d'));

set(handles.conf_slider,'Value',val);
update_figure(handles);

%temp = get(handles.conf_slider,'Value');
%set(hObject,'String',temp);
end

% --- Executes during object creation, after setting all properties.
function conf_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conf_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function conf_slider_Callback(hObject, eventdata, handles)
% hObject    handle to conf_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=round(hObject.Value);
hObject.Value=val;
set(handles.conf_number,'String',num2str(val,'%d'));
set(handles.A_electrode,'String',num2str(handles.A_elec(val),'%d'));
set(handles.M_electrode,'String',num2str(handles.M_elec(val),'%d'));
set(handles.N_electrode,'String',num2str(handles.N_elec(val),'%d'));
set(handles.B_electrode,'String',num2str(handles.B_elec(val),'%d'));

update_figure(handles);
end

function update_figure(handles)
A = str2double(get(handles.A_electrode,'String'));
M = str2double(get(handles.M_electrode,'String'));
N = str2double(get(handles.N_electrode,'String'));
B = str2double(get(handles.B_electrode,'String'));
axes(handles.axes1);
handles.current_ts=diff(handles.data3d_el_transient(A,[M N],:)-handles.data3d_el_transient(B,[M N],:));
if handles.switches.doPotentialsBase
    handles.current_base=diff(handles.basedata(A,[M N])-handles.basedata(B,[M N]));
end
if handles.switches.doPotentialsMoments
    m0 = diff(handles.data3d_el_moments(A,[M N],1)-handles.data3d_el_moments(B,[M N],1));
    m1 = diff(handles.data3d_el_moments(A,[M N],2)-handles.data3d_el_moments(B,[M N],2));
    handles.current_arrivaltime=m1/m0;
end
%timevector = datetime(datevec(seconds(handles.times)));
timevector = hours(handles.times/3600);
hold off
plot(timevector,reshape(handles.current_ts,[1,size(handles.current_ts,3)]),'kx-');%,'LineWidth',2);
if handles.switches.doPotentialsBase
    hold on
    plot([timevector(1),timevector(end)],[handles.current_base handles.current_base],'k--');%,'LineWidth',2);
end
if handles.switches.doPotentialsMoments
    hold on
    handles.current_arrivaltime=m1/m0;
    handles.current_arrivaltime = hours(handles.current_arrivaltime/3600);
    plot([handles.current_arrivaltime,handles.current_arrivaltime],[min(handles.current_ts) max(handles.current_ts)],'r-');
end

%xlim([0 handles.times]);
grid on;
ylabel('\Delta\phi in V');
xlabel('t (hh:mm)');
%ytickformat('%#3.2f')
xtickformat('hh:mm')
xticks(6*floor(timevector(1)/6,'hours'):hours(6):6*ceil(timevector(end)/6,'hours'))
set(handles.axes1,'XMinorTick','on','YMinorTick','on')
handles.axes1.XAxis.MinorTickValues = [6*floor(timevector(1)/6,'hours'):hours(1):6*ceil(timevector(end)/6,'hours')];
set(gca,'FontName','LM Roman 10','FontWeight','bold');%,'Color','none','FontSize',17)

axes(handles.axes5);
hold off
X = handles.X_elec;
Y = handles.Y_elec;
Z = handles.Z_elec;
scatter3(X,Y,Z,'o','MarkerEdgeColor',[0.4 0.4 0.4]);
grid minor
hold on
X = handles.X_elec([A]);
Y = handles.Y_elec([A]);
Z = handles.Z_elec([A]);
A = scatter3(X,Y,Z,'ko','filled','MarkerEdgeColor',[0 0 0]);
X = handles.X_elec([B]);
Y = handles.Y_elec([B]);
Z = handles.Z_elec([B]);
B = scatter3(X,Y,Z,'ro','filled','MarkerEdgeColor',[0 0 0]);
X = handles.X_elec([M]);
Y = handles.Y_elec([M]);
Z = handles.Z_elec([M]);
M = scatter3(X,Y,Z,'o','filled','MarkerFaceColor',[0 .9 .0],'MarkerEdgeColor',[0 0 0]);
X = handles.X_elec([N]);
Y = handles.Y_elec([N]);
Z = handles.Z_elec([N]);
N = scatter3(X,Y,Z,'o','filled','MarkerFaceColor',[0 .75 .9],'MarkerEdgeColor',[0 0 0]);
xlabel('x in m');
ylabel('y in m');
view([-16 72]);
legend([A,M,N,B],'A','M','N','B','Location','best');
set(gca,'FontName','LM Roman 10','FontWeight','bold');%,'Color','none','FontSize',13)
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key,'return') ||  strcmp(eventdata.Key,'tab')
    update_figure(handles);
    set(handles.conf_number,'String','manual')
end
end


% --- Executes on key press with focus on A_electrode and none of its controls.
function A_electrode_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function A_electrode_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function M_electrode_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function M_electrode_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function N_electrode_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function N_electrode_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function B_electrode_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function B_electrode_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function conf_number_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on key press with focus on A_electrode and none of its controls.
function conf_number_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to A_electrode (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on button press in set_A_electrode.
function set_A_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to set_A_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    axes(handles.axes5);
    dcm_obj = datacursormode(gcf);
    info_struct = getCursorInfo(dcm_obj);
    test_coord = info_struct.Position;
    
    index = find(handles.X_elec==test_coord(1) & handles.Y_elec==test_coord(2) & handles.Z_elec==test_coord(3));
    set(handles.A_electrode,'String',num2str(index));
end

update_figure(handles);
set(handles.conf_number,'String','manual')
end


% --- Executes on button press in set_M_electrode.
function set_M_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to set_M_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    axes(handles.axes5);
    dcm_obj = datacursormode(gcf);
    info_struct = getCursorInfo(dcm_obj);
    test_coord = info_struct.Position;
    index = find(handles.X_elec==test_coord(1) & handles.Y_elec==test_coord(2) & handles.Z_elec==test_coord(3));
    set(handles.M_electrode,'String',num2str(index));
end
update_figure(handles);
set(handles.conf_number,'String','manual')
end
% --- Executes on button press in set_N_electrode.
function set_N_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to set_N_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    axes(handles.axes5);
    dcm_obj = datacursormode(gcf);
    info_struct = getCursorInfo(dcm_obj);
    test_coord = info_struct.Position;
    index = find(handles.X_elec==test_coord(1) & handles.Y_elec==test_coord(2) & handles.Z_elec==test_coord(3));
    set(handles.N_electrode,'String',num2str(index));
end
update_figure(handles);
set(handles.conf_number,'String','manual')
end
% --- Executes on button press in set_B_electrode.
function set_B_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to set_B_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    axes(handles.axes5);
    dcm_obj = datacursormode(gcf);
    info_struct = getCursorInfo(dcm_obj);
    test_coord = info_struct.Position;
    index = find(handles.X_elec==test_coord(1) & handles.Y_elec==test_coord(2) & handles.Z_elec==test_coord(3));
    set(handles.B_electrode,'String',num2str(index));
end
update_figure(handles);
set(handles.conf_number,'String','manual')
end



function datafilename_Callback(hObject, eventdata, handles)
% hObject    handle to datafilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datafilename as text
%        str2double(get(hObject,'String')) returns contents of datafilename as a double
end

% --- Executes during object creation, after setting all properties.
function datafilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datafilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in printplots.
function printplots_Callback(hObject, eventdata, handles)
% hObject    handle to printplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sprintf('button');
filename1 = strcat('resistivity_A_',num2str(handles.A_electrode.String),'_M_',num2str(handles.M_electrode.String),'_N_',num2str(handles.N_electrode.String),'_B_',num2str(handles.B_electrode.String),'.pdf');
filename2 = strcat('configuration_A_',num2str(handles.A_electrode.String),'_M_',num2str(handles.M_electrode.String),'_N_',num2str(handles.N_electrode.String),'_B_',num2str(handles.B_electrode.String),'.pdf');
prompt = {'Filename Development','Filename Configuration'};
title = 'Plot';
dims = [1 35];
definput = {filename1,filename2};
answer = inputdlg(prompt,title,dims,definput);

figure_handle = isolate_axes(handles.axes1);
handles.axes1.SortMethod='ChildOrder';
export_fig(figure_handle, answer{1},'-p10','-transparent','-linecaps','-painters'); % Add the correct print options here
close(figure_handle);

figure_handle = isolate_axes(handles.axes5);
handles.axes5.SortMethod='ChildOrder';
export_fig(figure_handle, answer{2},'-p10','-transparent','-linecaps','-painters'); % Add the correct print options here
close(figure_handle);
end
