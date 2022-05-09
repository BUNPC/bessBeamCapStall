function varargout = capStall(varargin)
% CAPSTALL MATLAB code for capStall.fig
%      CAPSTALL, by itself, creates a new CAPSTALL or raises the existing
%      singleton*.ans
%
%      H = CAPSTALL returns the handle to a new CAPSTALL or the handle to
%      the existing singleton*.
%
%      CAPSTALL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAPSTALL.M with the given input arguments.
%
%      CAPSTALL('Property','Value',...) creates a new CAPSTALL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before capStall_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to capStall_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu_loaddata.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help capStall

% Last Modified by GUIDE v2.5 09-May-2022 00:56:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @capStall_OpeningFcn, ...
                   'gui_OutputFcn',  @capStall_OutputFcn, ...
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


% --- Executes just before capStall is made visible.
function capStall_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to capStall (see VARARGIN)

% Choose default command line output for capStall
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes capStall wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = capStall_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_loadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

% clear Data strutre
if exist('Data','var')
    if isstruct(Data)
        names = fieldnames(Data);
        Data = rmfield(Data,names);
    end
end

[filename,pathname] = uigetfile({'*.tiff;*.tif'},'Please select the Angiogram Data');
if filename == 0
    return
end
handles.text_Imagepath.String = [pathname filename];
handles.text_Imagepath.Visible = 'on';
h = waitbar(0,'Please wait... loading the data');
[~,~,ext] = fileparts(filename);

if strcmp(ext,'.mat')
    temp = load([pathname filename]);
    fn = fieldnames(temp);
    angio = temp.(fn{1});
elseif  strcmp(ext,'.tiff') || strcmp(ext,'.tif')
    info = imfinfo([pathname filename]);
    for u = 1:length(info)
        if u == 1
            temp = imread([pathname filename],1);
            angio = zeros([length(info) size(temp)]);
            angio(u,:,:) = temp;
        else
            angio(u,:,:) = imread([pathname filename],u);
        end
    end
end
close(h)


Data.pathname = pathname;
Data.I = permute(angio,[2 3 1]);
Data.Volume = mean(angio,1);
Vz = size(Data.Volume,1);
[Sx, Sy, Sz] = size(Data.I);

set(handles.edit_Xmin,'String',num2str(1));
set(handles.edit_Ymin,'String',num2str(1));
set(handles.edit_Xmax,'String',num2str(Sx));
set(handles.edit_Ymax,'String',num2str(Sy));
set(handles.text_Xrange,'String',['Xrange(1-' num2str(Sx) ')']);
set(handles.text_Yrange,'String',['Yrange(1-' num2str(Sy) ')']);

set(handles.slider_movedata,'max',Sz);
set(handles.slider_movedata,'min',1);
set(handles.slider_movedata,'Value',1);
set(handles.slider_movedata,'SliderStep',[1/(Sz-1), 10/(Sz-1)]);

img = Data.I(:,:,1);
set(handles.edit_MaxI,'String',num2str(max(img(:))));
set(handles.edit_MinI,'String',num2str(min(img(:))));

draw(hObject, eventdata, handles);

function draw(hObject, eventdata, handles)

global Data
I = Data.I;
ii = str2double(get(handles.edit_volnumber,'string'));
MinI = str2double(get(handles.edit_MinI,'String'));
MaxI = str2double(get(handles.edit_MaxI,'String'));
Xmin = max(str2double(get(handles.edit_Xmin,'String')),1);
Ymin = max(str2double(get(handles.edit_Ymin,'String')),1);
Xmax = str2double(get(handles.edit_Xmax,'String'));
Ymax = str2double(get(handles.edit_Ymax,'String'));

if isfield(Data,'Int_ts')
    cseg = str2double(get(handles.edit_segno,'string'));
else 
    cseg = 0;
end

axes(handles.axes1)
colormap('gray');
h = imagesc(I(:,:,ii),[MinI MaxI]);
axis image;
xlim([Xmin Xmax])
ylim([Ymin Ymax])
hold on
if isfield(Data,'Cap')
    for u = 1:size(Data.Cap,1)
        if Data.Cap(u,2) >= Xmin & Data.Cap(u,2) <= Xmax & Data.Cap(u,1) >= Ymin & Data.Cap(u,1) <= Ymax
            if u == cseg
                hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','r','FontSize',10);
            else
                hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','g','FontSize',10);
            end
            set(hpt,'ButtonDownFcn', sprintf('Cap_Stall_deletept(%d)',u) );
        end
    end
end
hold off
set(h, 'ButtonDownFcn', {@axes1_ButtonDown, handles},'BusyAction','cancel');
set(gcf, 'WindowScrollWheelFcn', {@axes_WindowScrollWheelFcn, handles},'Interruptible','off','BusyAction','cancel');
% set(gcf, 'WindowKeyPressFcn', {@figure_WindowKeyPressFcn, handles},'Interruptible','off','BusyAction','cancel');

seg_no = str2double(get(handles.edit_segno,'string'));
frame_no = str2double(get(handles.edit_volnumber,'string'));
msg_display = '';
if strcmpi(get(handles.menu_validateStalls,'Checked'), 'on')
    if isfield(Data,'seg')
        if isfield(Data,'StallingMatrix')
            msg_display = 'Human:';
            set(handles.text_stallinfoName,'String',msg_display, 'fontsize', 10)
            if Data.StallingMatrix(seg_no,frame_no) == 1
                msg_display = 'Stall';
                color = 'r';
            else
                msg_display = 'Not a Stall';
                color = 'k';
            end
            set(handles.text_stallinfo,'String',msg_display, 'fontsize', 10)
            set(handles.text_stallinfo,'ForegroundColor',color)
        end

        if isfield(Data,'AutoStallingMatrix')
            msg_display = 'Auto:';
            set(handles.text_AutostallinfoName,'String',msg_display, 'fontsize', 10)
            if Data.AutoStallingMatrix(seg_no,frame_no) == 1
                msg_display = 'Stall';
                color = 'g';
            else
                msg_display = 'Not a Stall';
                color = 'k';
            end
            set(handles.text_Autostallinfo,'String',msg_display, 'fontsize', 10)
            set(handles.text_Autostallinfo,'ForegroundColor',color)
        end

        if isfield(Data,'GTStallingMatrix')
            msg_display = 'GT:';
            set(handles.text_GTstallinfoName,'String',msg_display, 'fontsize', 10)
            if Data.GTStallingMatrix(seg_no,frame_no) == 1
                msg_display = 'Stall';
                color = 'b';
            elseif Data.GTStallingMatrix(seg_no,frame_no) == 2
                msg_display = 'Questionable';
                color = 'm';
            elseif Data.GTStallingMatrix(seg_no,frame_no) == 3
                msg_display = 'Skipped';
                color = 'c';
            else
                msg_display = 'Not a Stall';
                color = 'k';
            end
            set(handles.text_GTstallinfo,'String',msg_display, 'fontsize', 10)
            set(handles.text_GTstallinfo,'ForegroundColor',color)
        end
    end
end

if handles.checkbox_Axe1displayCurrentSegment.Value == 1
    if isfield(Data.seg,'pos')
        hold on
        if isfield(Data.seg(seg_no),'frame_seg_pos')
            plot(Data.seg(seg_no).frame_seg_pos(:,2,frame_no),Data.seg(seg_no).frame_seg_pos(:,1,frame_no),'r.','markersize',16);
        else
            plot(Data.seg(seg_no).pos(:,1),Data.seg(seg_no).pos(:,2),'r.','markersize',16);
        end
        hold off
    else
        handles.checkbox_Axe1displayCurrentSegment.Value = 0;
        error('Warning: No segments available')
    end
end


axes(handles.axes2)
colormap('gray');
if isfield(Data,'seg')
    delete(handles.axes2.Children)
    if isfield(Data.seg,'LRimage') || isfield(Data.seg,'pos')
        hold on
        [Sx, Sy, ~] = size(Data.I);
        seg_pos = Data.seg(seg_no).pos;
        deltaX = str2double(get(handles.edit_deltaX,'String'));
        deltaY = str2double(get(handles.edit_deltaY,'String'));
        Xmin = min(seg_pos(:,1))-deltaX;
        Xmax = max(seg_pos(:,1))+deltaX;
        Ymin = min(seg_pos(:,2))-deltaY;
        Ymax = max(seg_pos(:,2))+deltaY;
        Xmin = min(max(Xmin,1),Sx);
        Xmax = min(max(Xmax,Xmin+1),Sx);
        Ymin = min(max(Ymin,1),Sy);
        Ymax = min(max(Ymax,Ymin+1),Sy);
        h2 = imagesc(I(:,:,ii),[MinI MaxI]);
        axis image;
        xlim([Xmin Xmax])
        ylim([Ymin Ymax])
        if handles.checkbox_Axe2displayCurrentSegment.Value == 1
            if isfield(Data.seg,'pos')
                hold on
                if isfield(Data.seg(seg_no),'frame_seg_pos')
                    plot(Data.seg(seg_no).frame_seg_pos(:,2,frame_no),Data.seg(seg_no).frame_seg_pos(:,1,frame_no),'r.','markersize',3,'Marker','o');
                else
                    plot(Data.seg(seg_no).pos(:,1),Data.seg(seg_no).pos(:,2),'r.','markersize',16);
                end
                hold off
            end
        end
        hold off
    else
        colormap('gray');
        h2 = imagesc(I(:,:,ii),[MinI MaxI]);
        axis image;
    end
else
    colormap('gray');
    h2 = imagesc(I(:,:,ii),[MinI MaxI]);
    axis image;
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
end

if isfield(Data,'seg') && isfield(Data.seg,'LRimage')
    
    jj = str2double(get(handles.edit_segno,'string'));
    axes(handles.axes3)
    hp = imagesc(Data.seg(jj).LRimage');
    axis image
    hold on
    x = xlim;
    handles.axes3.YAxis.Limits = [0,351];
    handles.axes3.XTick = [];
    handles.axes3.YTick = linspace(0,350,36);
    y = [ii ii];
    line(x,y, 'LineWidth', 3, 'Color','r');
    if isfield(Data,'StallingMatrix')
        xidx = find(Data.StallingMatrix(jj,:) == 1);
        yidx = size(Data.seg(jj).LRimage,1)/2*ones(1,length(xidx));
        text(yidx,xidx,'*','Color','red','FontSize',10);
    end
    if isfield(Data,'AutoStallingMatrix')
        xidx = find(Data.AutoStallingMatrix(jj,:) == 1);
        yidx = size(Data.seg(jj).LRimage,1)/3*ones(1,length(xidx));
        text(yidx,xidx,'*','Color','green','FontSize',10);
    end
    
    if isfield(Data,'GTStallingMatrix')
        xidx = find(Data.GTStallingMatrix(jj,:) == 1);
        yidx = size(Data.seg(jj).LRimage,1)/1.5*ones(1,length(xidx));
        text(yidx,xidx,'*','Color','blue','FontSize',10);
        
        xidx = find(Data.GTStallingMatrix(jj,:) == 2);
        yidx = size(Data.seg(jj).LRimage,1)/1.5*ones(1,length(xidx));
        text(yidx,xidx,'*','Color','magenta','FontSize',10);
        
        xidx = find(Data.GTStallingMatrix(jj,:) == 3);
        yidx = size(Data.seg(jj).LRimage,1)/1.5*ones(1,length(xidx));
        text(yidx,xidx,'*','Color','cyan','FontSize',10);
    end
    hold off
    handles.text_SegLengthDisplay.String = [ num2str(x(2)-x(1)) ' Pixels' ];
    set(handles.axes3, 'ButtonDownFcn', {@axes3_ButtonDownFcn, handles});
end

Data.handles =  handles;
Data.hObject =  hObject;
Data.eventdata = eventdata;


function axes3_ButtonDownFcn(hObject, eventdata, handles)

global Data

parent = (get(hObject, 'Parent'));
mouseclick = get(parent, 'SelectionType');
if strcmp(mouseclick,'normal')
    ii = round(eventdata.IntersectionPoint(2));
    %ii = min(max(ii,1),size(Data.I,3));
    if ii <= 350 && ii >=1
        set(handles.edit_volnumber,'string',num2str(ii));
        set(handles.slider_movedata,'Value',ii);
        draw(hObject, eventdata, handles);
        if isfield(Data,'sliderobject')
            uicontrol(Data.sliderobject);
        end
    end

end
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobject')
    uicontrol(Data.sliderobject);
end

function axes_WindowScrollWheelFcn(hObject, eventdata, handles)

global Data

[I_x,I_y,~] = size(Data.I);
[~,V_x,V_y] = size(Data.Volume);
axis1pos = get(handles.axes1, 'CurrentPoint');
axis1pos = axis1pos(1,1:2);
axis2pos = get(handles.axes2,'CurrentPoint');
axis2pos = axis2pos(1,1:2);
axis3pos = get(handles.axes3,'CurrentPoint');
axis3pos = axis3pos(1,1:2);
a3x = handles.axes3.XLim;
a3y = handles.axes3.YLim;
if isfield(Data,'Int_ts')
    cseg = str2double(get(handles.edit_segno,'string'));
else 
    cseg = 0;
end
% if  (axis3pos(1) >=a3x(1) && axis3pos(1) <= a3x(2) && axis3pos(2) >=a3y(1) && axis3pos(2) <=a3y(2))
    if eventdata.VerticalScrollCount > 0
        ii = str2num(get(handles.edit_volnumber,'string'));
        ii = min(max(ii+1,1),size(Data.I,3));
        set(handles.edit_volnumber,'string',num2str(ii));
        set(handles.slider_movedata,'Value',ii);
    end
    if eventdata.VerticalScrollCount < 0
        ii = str2num(get(handles.edit_volnumber,'string'));
        ii = min(max(ii-1,1),size(Data.I,3));
        set(handles.edit_volnumber,'string',num2str(ii));
        set(handles.slider_movedata,'Value',ii);
    end
    draw(hObject, eventdata, handles);

function edit_volnumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_volnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_volnumber as text
%        str2double(get(hObject,'String')) returns contents of edit_volnumber as a double

global Data
ii = str2num(get(handles.edit_volnumber,'string'));
ii = min(max(ii,1),size(Data.I,3));
set(handles.edit_volnumber,'string',num2str(ii));
set(handles.slider_movedata,'Value',ii);
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobject')
    uicontrol(Data.sliderobject);
end

% --- Executes during object creation, after setting all properties.
function edit_volnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_volnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_moveleft.
function pushbutton_moveleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2num(get(handles.edit_volnumber,'string'));
ii = min(max(ii-1,1),size(Data.I,3));
set(handles.edit_volnumber,'string',num2str(ii));
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobject')
    uicontrol(Data.sliderobject);
end

% --- Executes on button press in pushbutton_moveright.
function pushbutton_moveright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_moveright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2num(get(handles.edit_volnumber,'string'));
ii = min(max(ii+1,1),size(Data.I,3));
set(handles.edit_volnumber,'string',num2str(ii));
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobject')
    uicontrol(Data.sliderobject);
end


% --- Executes on key press with focus on edit_volnumber and none of its controls.
function edit_volnumber_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_volnumber (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

disp('here');


% --- Executes on slider movement.
function slider_movedata_Callback(hObject, eventdata, handles)
% hObject    handle to slider_movedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global Data

ii = round(get(hObject,'Value'));
ii = min(max(ii+1,1),size(Data.I,3));
set(handles.edit_volnumber,'string',num2str(ii)); 
Data.sliderobject = hObject;
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider_movedata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_movedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
global Data
Data.sliderobject = hObject;


function axes1_ButtonDown(hObject, eventdata, handles)

global Data

parent = (get(hObject, 'Parent'));
pts = round(get(parent, 'CurrentPoint'));
y = pts(1,1);
x = pts(1,2);
if isfield(Data,'Cap')
    Data.Cap = [Data.Cap;x y];
else
    Data.Cap = [x y];
end

draw(hObject, eventdata, handles);
if isfield(Data,'sliderobject')
    uicontrol(Data.sliderobject);
end

% --- Executes on slider movement.
function slider_moveinZ_Callback(hObject, eventdata, handles)
% hObject    handle to slider_moveinZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global Data

Vz = size(Data.Volume,1);
ii = Vz-round(get(handles.slider_moveinZ,'Value'))+1;
jj = str2double(get(handles.edit_MIPnofframes,'String'));
ii = min(max(ii,1),Vz);
jj = min(max(jj,1),Vz-ii+1);
set(handles.edit_MIPstartidx,'String',num2str(ii));
set(handles.edit_MIPnofframes,'String',num2str(jj));
Data.sliderobjectZ = hObject;
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider_moveinZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_moveinZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
global Data
Data.sliderobjectZ = hObject;



function edit_MIPstartidx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MIPstartidx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MIPstartidx as text
%        str2double(get(hObject,'String')) returns contents of edit_MIPstartidx as a double

global Data

Vz = size(Data.Volume,1);
ii = str2double(get(handles.edit_MIPstartidx,'String'));
jj = str2double(get(handles.edit_MIPnofframes,'String'));
ii = min(max(ii,1),Vz);
jj = min(max(jj,1),Vz-ii+1);
set(handles.edit_MIPstartidx,'String',num2str(ii));
set(handles.edit_MIPnofframes,'String',num2str(jj));
set(handles.slider_moveinZ,'Value',Vz-ii+1);
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobjectZ')
    uicontrol(Data.sliderobjectZ);
end


% --- Executes during object creation, after setting all properties.
function edit_MIPstartidx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MIPstartidx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MIPnofframes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MIPnofframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MIPnofframes as text
%        str2double(get(hObject,'String')) returns contents of edit_MIPnofframes as a double

global Data

Vz = size(Data.Volume,1);
ii = str2double(get(handles.edit_MIPstartidx,'String'));
jj = str2double(get(handles.edit_MIPnofframes,'String'));
ii = min(max(ii,1),Vz);
jj = min(max(jj,1),Vz-ii+1);
set(handles.edit_MIPstartidx,'String',num2str(ii));
set(handles.edit_MIPnofframes,'String',num2str(jj));
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobjectZ')
    uicontrol(Data.sliderobjectZ);
end



% --- Executes during object creation, after setting all properties.
function edit_MIPnofframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MIPnofframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_SegAnalysis.
function pushbutton_SegAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SegAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

img = squeeze(Data.Volume);
img_VE = fibermetric(img,6);
img_seg = imbinarize(img_VE,'adaptive');

Data.Volume_seg = img_seg;

angio_2 = zeros([size(Data.Volume_seg) 2]);
angio_2(:,:,1) = Data.Volume_seg;
vessel_mask = logical(angio_2);

vessel_skl = bwskel(vessel_mask);

vessel_graph = fun_skeleton_to_graph(vessel_skl);
vessel_mask_dt = bwdist(~vessel_mask);


% Size of the angiogram. It will help to convert indeces to subscripts
angio_size = size(vessel_mask);

% find length of nodes to allocate size
nodes_ind_count = length(vessel_graph.node.cc_ind)+length(vessel_graph.link.pos_ind);
nodes_ind = zeros(1,nodes_ind_count);

% find length of edges to allocate size
edges_ind_count = 0;
% link_cc_ind = zeros(size(vessel_graph.link.cc_ind));
for u = 1:length(vessel_graph.node.connected_link_label)
    edges_ind_count = edges_ind_count+length(vessel_graph.node.connected_link_label{u});
end
for  u = 1:length(vessel_graph.link.cc_ind)
    edges_ind_count = edges_ind_count+length(vessel_graph.link.cc_ind{u})-1;
end
% edges_ind_count = edges_ind_count+length(vessel_graph.link.pos_ind);
% idx = find(link_cc_ind == 0);
% for u = 1:length(idx)
%     
% end
edges_ind = zeros(edges_ind_count,2);

% assign nodes and edges
node_idx = 1;
edge_idx = 1;
tttt = [];
link_cc_ind = zeros(size(vessel_graph.link.cc_ind));
for u = 1:length(vessel_graph.node.cc_ind)
    nodes_ind(node_idx) = vessel_graph.node.cc_ind{u}(1);
    tttt = [tttt node_idx];
    temp_node = nodes_ind(node_idx);
    for v = 1:length(vessel_graph.node.connected_link_label{u})
        connected_link = vessel_graph.node.connected_link_label{u}(v);
        connected_link_endnodes = [vessel_graph.link.cc_ind{connected_link}(1) vessel_graph.link.cc_ind{connected_link}(end)];
        [n1,n2,n3] = ind2sub(angio_size,temp_node);
        for w = 1:2
            [l1,l2,l3] = ind2sub(angio_size,connected_link_endnodes(w));
            d(w) = sqrt((n1-l1)^2+(n2-l2)^2+(n3-l3)^2);
        end
        [~,min_idx] = min(d);
        edges_ind(edge_idx,:) = [temp_node connected_link_endnodes(min_idx(1))];
        if link_cc_ind(connected_link) == 0
            link_length = length(vessel_graph.link.cc_ind{connected_link});
            edges_ind(edge_idx+1:edge_idx+link_length-1,1) = vessel_graph.link.cc_ind{connected_link}(1:end-1);
            edges_ind(edge_idx+1:edge_idx+link_length-1,2) = vessel_graph.link.cc_ind{connected_link}(2:end);
            edge_idx = edge_idx+link_length;
            nodes_ind(node_idx+1:node_idx+link_length) = vessel_graph.link.cc_ind{connected_link};
            isa = ismember(nodes_ind(node_idx+1:node_idx+link_length),0);
            tttt = [tttt node_idx+1:node_idx+link_length];
            node_idx = node_idx+link_length;
            link_cc_ind(connected_link) = 1;
        else
            edge_idx = edge_idx+1;
        end
    end
    node_idx = node_idx+1;
end

idx = find(link_cc_ind == 0);
for u = 1:length(idx)
    link_length = length(vessel_graph.link.cc_ind{idx(u)})
    edges_ind(edge_idx+1:edge_idx+link_length-1,1) = vessel_graph.link.cc_ind{idx(u)}(1:end-1);
    edges_ind(edge_idx+1:edge_idx+link_length-1,2) = vessel_graph.link.cc_ind{idx(u)}(2:end);
    edge_idx = edge_idx+link_length;
    nodes_ind(node_idx+1:node_idx+link_length) = vessel_graph.link.cc_ind{idx(u)};
    isa = ismember(nodes_ind(node_idx+1:node_idx+link_length),0);
    tttt = [tttt node_idx+1:node_idx+link_length];
    node_idx = node_idx+link_length;
end

% nodes = zeros(length(nodes_ind),3);
[n1 n2 n3] = ind2sub(angio_size,nodes_ind);
nodes =[n1' n2' n3'];
edges = zeros(size(edges_ind));

for u = 1:size(edges_ind,1)
    u
    edges(u,1) = find(nodes_ind == edges_ind(u,1));
    edges(u,2) = find(nodes_ind == edges_ind(u,2));
end

Graph.nodes = nodes;
Graph.edges = edges;

sameEdgeIdx = [];
for u = 1:size(Graph.edges,1)
    if Graph.edges(u,1) == Graph.edges(u,2)
        sameEdgeIdx = [sameEdgeIdx; u];
    end
end
Graph.edges(sameEdgeIdx,:) = [];

temp = Graph.nodes(:,2);
Graph.nodes(:,2) = Graph.nodes(:,1);
Graph.nodes(:,1) = temp;
Data.Graph = Graph;

Data.Graph.segInfo = nodeGrps_vesSegment(Data.Graph.nodes, Data.Graph.edges);
% n_frames = size(Data.I,3); 
[bX,bY,n_frames] = size(Data.I);
for vv = 1:size(Data.Cap,1)
    vv
    dist = sqrt(sum((Data.Graph.nodes(:,1:2)-[Data.Cap(vv,2) Data.Cap(vv,1)]).^2,2));
    [min_dist,min_idx] = min(dist);
    seg_no = Data.Graph.segInfo.nodeSegN(min_idx);
    seg_nodes = find(Data.Graph.segInfo.nodeSegN == seg_no);
    seg(vv).pos = Data.Graph.nodes(seg_nodes,:);
    seg(vv).mask = Data.Graph.nodes(seg_nodes,:);
    
    min_seg_x = min(seg(vv).pos(:,2));
    max_seg_x = max(seg(vv).pos(:,2));
    min_seg_y = min(seg(vv).pos(:,1));
    max_seg_y = max(seg(vv).pos(:,1));
    min_x = min(max(min_seg_x-7,1),bY);
    max_x = min(max(max_seg_x+7,1),bY);
    min_y = min(max(min_seg_y-7,1),bX);
    max_y = min(max(max_seg_y+7,1),bX);
    mean_cropped_image = squeeze(Data.Volume(1,min_x:max_x,min_y:max_y));
    if max_x-min_x+1 < 16
        mean_cropped_image = [mean_cropped_image; zeros(15-(max_x-min_x), size(mean_cropped_image,2))];
    end
    if max_y-min_y+1 < 16
        mean_cropped_image = [mean_cropped_image zeros(size(mean_cropped_image,1), 15-(max_y-min_y))];
    end
    mean_seg_pos = [seg(vv).pos(:,2)-min_x seg(vv).pos(:,1)-min_y seg(vv).pos(:,3)];
    mean_mask = zeros(size(mean_cropped_image));
    [optimizer, metric] = imregconfig('monomodal');
    LRimage = zeros(size(seg(vv).pos,1),n_frames);
    frame_seg_pos = zeros(size(seg(vv).pos,1),2,n_frames);
    for ww =1:n_frames
%         ww
        frame_cropped_image = Data.I(min_x:max_x,min_y:max_y,ww);
        current_frame = squeeze(Data.I(:,:,ww));
        if max_x-min_x+1 < 16
            frame_cropped_image = [frame_cropped_image; zeros(15-(max_x-min_x), size(frame_cropped_image,2))];
        end
        if max_y-min_y+1 < 16
            frame_cropped_image = [frame_cropped_image zeros(size(frame_cropped_image,1), 15-(max_y-min_y))];
        end
        lastwarn('')
        tform = imregtform(mean_cropped_image', frame_cropped_image', 'rigid', optimizer, metric);
        mean_seg_pos = [seg(vv).pos(:,2)-min_x seg(vv).pos(:,1)-min_y seg(vv).pos(:,3)];
        [warnMsg, warnId] = lastwarn;
        if isempty(warnMsg)
            [points_x, points_y] = transformPointsForward(tform,mean_seg_pos(:,1),mean_seg_pos(:,2));
        else
            points_x = mean_seg_pos(:,1); 
            points_y = mean_seg_pos(:,2);
        end
        frame_seg_pos_x =points_x+min_x ;
        frame_seg_pos_y = points_y+min_y;
        frame_seg_pos_x(frame_seg_pos_x<1) = 1;
        frame_seg_pos_x(frame_seg_pos_x>bY) = bY;
        frame_seg_pos_y(frame_seg_pos_y<1) = 1;
        frame_seg_pos_y(frame_seg_pos_y>bX) = bX;
%         ind = round(sub2ind(size(current_frame),frame_seg_pos_x',frame_seg_pos_y'));
%         LRimage(:,ww) = current_frame(ind);
        frame_seg_pos(:,:,ww) = [frame_seg_pos_x frame_seg_pos_y];
        for ll = 1:size(seg(vv).pos,1)
            LRimage(ll,ww) = current_frame(round(frame_seg_pos_x(ll)),round(frame_seg_pos_y(ll)));
        end
%         moving_reg = imwarp(mean_cropped_image, tform, 'OutputView',imref2d(size(cropped_image)));
    end
    seg(vv).LRimage = LRimage;
    seg(vv).frame_seg_pos = frame_seg_pos;
end
Data.seg = seg;
Int_ts = zeros(length(seg),n_frames);
for kk = 1:length(seg)
    Int_ts(kk,:) = mean(seg(kk).LRimage,1);
end

Data.Int_ts = Int_ts;

MinI = str2double(get(handles.edit_MinI,'String'));
MaxI = str2double(get(handles.edit_MaxI,'String'));
figure; colormap('gray');
imagesc(squeeze(max(Data.Volume,[],1)),[MinI MaxI]); 
hold on
for uu = 1:length(seg)
    plot(seg(uu).pos(:,1),seg(uu).pos(:,2),'r.','markersize',12);
end
hold off
axis image

% StallingMatrix for manual annotation
if isfield(Data,'StallingMatrix')
   Int_rows = size(Data.Int_ts,1);
   Stall_rows = size(Data.StallingMatrix,1);
   stall_cols = size(Data.StallingMatrix,2);
   if Int_rows > Stall_rows
       Data.StallingMatrix = [Data.StallingMatrix; zeros([Int_rows-Stall_rows stall_cols])];
   end
else
   Data.StallingMatrix = zeros(length(seg),n_frames);
end

% Find stalls using cross correlation
smoothFactor = 1;
corrThresh = 0.85;
Data.AutoStallingMatrix = zeros(length(seg),n_frames);
for i = 1:length(seg)
    crossVals = correlateLT(seg(i).LRimage,smoothFactor);
    idx = crossVals>=corrThresh;
    Data.AutoStallingMatrix(i,idx) = 1;
end

% create or update GT stalling matrix
if isfield(Data,'GTStallingMatrix')
    Int_rows = length(Data.seg);
    Stall_rows = size(Data.GTStallingMatrix,1);
    stall_cols = size(Data.GTStallingMatrix,2);
    if Int_rows > Stall_rows
       Data.GTStallingMatrix = [Data.GTStallingMatrix; zeros([Int_rows-Stall_rows stall_cols])];
    end
else
    Data.GTStallingMatrix = zeros(length(seg),n_frames);
    idx = find(Data.StallingMatrix == 1 & Data.AutoStallingMatrix==1);
    Data.GTStallingMatrix(idx) = 1;
end

% create or update Validation Flag
if isfield(Data,'ValidationFlag')
    Int_rows = length(Data.seg);
    Stall_rows = size(Data.ValidationFlag,1);
    stall_cols = size(Data.ValidationFlag,2);
    if Int_rows > Stall_rows
        update = Data.StallingMatrix(Stall_rows+1:end,:) == Data.AutoStallingMatrix(Stall_rows+1:end,:);
        Data.ValidationFlag = [Data.ValidationFlag; update];
    end
else
    Data.ValidationFlag = zeros(length(Data.seg),n_frames);
    idx = find(Data.StallingMatrix == 1 & Data.AutoStallingMatrix==1);
    Data.ValidationFlag(idx) = 1;
    idx = find(Data.StallingMatrix == 0 & Data.AutoStallingMatrix==0);
    Data.ValidationFlag(idx) = 1;
end
draw(hObject, eventdata, handles);
makeSegLengthHistogram(handles);


function crossCorrVals = correlateLT(LRimage,smoothFactor)
    LRimage_mean = mean(LRimage ,2);
    LRimage_contrast = LRimage-LRimage_mean;
    maxLag =2;
    
    cc = [];
    for u = 1:size(LRimage_contrast,2)-1
        R = xcorr(LRimage_contrast(:,u),LRimage_contrast(:,u+1),maxLag,'normalized');
        cc = [cc ; max(R)];
    end
    f_coeff = ones(1,smoothFactor)/smoothFactor;
    crossCorrVals = filter(f_coeff, 1, cc);
   

function edit_segno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_segno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_segno as text
%        str2double(get(hObject,'String')) returns contents of edit_segno as a double

global Data
jj = str2double(get(handles.edit_segno,'string'));
jj = min(max(jj,1),length(Data.seg));
set(handles.edit_segno,'string',num2str(jj));
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_segno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_segno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_prevseg.
function pushbutton_prevseg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_prevseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
jj = str2double(get(handles.edit_segno,'string'));
jj = min(max(jj-1,1),length(Data.seg));
set(handles.edit_segno,'string',num2str(jj));
draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_nextseg.
function pushbutton_nextseg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nextseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
jj = str2double(get(handles.edit_segno,'string'));
jj = min(max(jj+1,1),length(Data.seg));
set(handles.edit_segno,'string',num2str(jj));
draw(hObject, eventdata, handles);

% --- Executes on button press in radiobutton_showseg.
function radiobutton_showseg_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_showseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_showseg

draw(hObject, eventdata, handles);



function edit_MaxI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxI as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxI as a double

draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_MaxI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MinI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinI as text
%        str2double(get(hObject,'String')) returns contents of edit_MinI as a double

draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_MinI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_loadData_Callback(hObject, eventdata, handles)
% hObject    handle to menu_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data;
handles.checkbox_Axe1displayCurrentSegment.Value = 0;
delete(handles.axes3.Children)
delete(handles.axes4.Children)

[filename, pathname] = uigetfile;
if filename == 0
    return
end
temp_struct = load([pathname filename]);
handles.text_Datapath.String = [pathname filename];
handles.text_Datapath.Visible = 'on';
if isfield(temp_struct,'Cap')
    Data.Cap = temp_struct.Cap;
end
if isfield(temp_struct,'pts')
    Data.pts = temp_struct.pts;
end
if isfield(temp_struct,'seg')
    Data.seg = temp_struct.seg;
    if isfield(Data.seg,'LRimage')
        makeSegLengthHistogram(handles);
    end
end
if isfield(temp_struct,'Int_ts')
    Data.Int_ts = temp_struct.Int_ts;
end
if isfield(temp_struct,'StallingMatrix')
   Data.StallingMatrix = temp_struct.StallingMatrix;
end
if isfield(temp_struct,'AutoStallingMatrix')
   Data.AutoStallingMatrix = temp_struct.AutoStallingMatrix;
end 
if isfield(temp_struct,'GTStallingMatrix')
   Data.GTStallingMatrix = temp_struct.GTStallingMatrix; 
   handles.checkbox_Axe1displayCurrentSegment.Value = 1;
end 
if isfield(temp_struct,'ValidationFlag')
   Data.ValidationFlag = temp_struct.ValidationFlag; 
end 
if isfield(temp_struct,'segAnalysis')
    Data.segAnalysis = temp_struct.segAnalysis; 
end
draw(hObject, eventdata, handles);


% --------------------------------------------------------------------
function menu_saveresults_Callback(hObject, eventdata, handles)
% hObject    handle to menu_saveresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

[filename, pathname] = uiputfile('*.mat');

if isfield(Data,'Cap')
    Cap = Data.Cap;
    save([pathname filename],'Cap');
end
if isfield(Data,'pts')
    pts = Data.pts;
    save([pathname filename],'pts','-append');
end
if isfield(Data,'seg')
    seg = Data.seg;
    save([pathname filename],'seg','-append');
end
if isfield(Data,'Int_ts')
    Int_ts=  Data.Int_ts;
    save([pathname filename],'Int_ts','-append');
end
if isfield(Data,'StallingMatrix')
    StallingMatrix = Data.StallingMatrix;
    save([pathname filename],'StallingMatrix','-append');
end
if isfield(Data,'AutoStallingMatrix')
    AutoStallingMatrix = Data.AutoStallingMatrix;
    save([pathname filename],'AutoStallingMatrix','-append');
end 
if isfield(Data,'GTStallingMatrix')
    GTStallingMatrix = Data.GTStallingMatrix;
    save([pathname filename],'GTStallingMatrix','-append');
end 
if isfield(Data,'ValidationFlag')
    ValidationFlag = Data.ValidationFlag;
    save([pathname filename],'ValidationFlag','-append');
end 
if isfield(Data,'segAnalysis')
    segAnalysis = Data.segAnalysis;
    save([pathname filename],'segAnalysis','-append');
end



function edit_Xmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xmin as text
%        str2double(get(hObject,'String')) returns contents of edit_Xmin as a double

global Data

[Sx,~, ~] = size(Data.I);
Xmin = str2double(get(hObject,'String'));
Xmin = min(max(Xmin,1),Sx);
set(handles.edit_Xmin,'String',num2str(Xmin))
draw(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Xmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xmax as text
%        str2double(get(hObject,'String')) returns contents of edit_Xmax as a double
global Data

[Sx,~, ~] = size(Data.I);
Xmin = str2double(get(handles.edit_Xmin,'String'));
Xmax = str2double(get(hObject,'String'));
Xmax = min(max(Xmax,Xmin+1),Sx);
set(handles.edit_Xmax,'String',num2str(Xmax))
draw(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ymin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ymin as text
%        str2double(get(hObject,'String')) returns contents of edit_Ymin as a double

global Data

[~,Sy, ~] = size(Data.I);
Ymin = str2double(get(hObject,'String'));
Ymin = min(max(Ymin,1),Sy);
set(handles.edit_Ymin,'String',num2str(Ymin))

draw(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ymax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ymax as text
%        str2double(get(hObject,'String')) returns contents of edit_Ymax as a double

global Data

[~,Sy, ~] = size(Data.I);
Ymin = str2double(get(handles.edit_Ymin,'String'));
Ymax = str2double(get(hObject,'String'));
Ymax = min(max(Ymax,Ymin+1),Sy);
set(handles.edit_Ymax,'String',num2str(Ymax))

draw(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_selectZoomInArea.
function pushbutton_selectZoomInArea_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectZoomInArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data

[Sx, Sy, ~] = size(Data.I);
k = waitforbuttonpress;
if k==0
    point1 = get(handles.axes1,'CurrentPoint');     % button down detected
    finalRect = rbbox;                              % return figure units
    point2 = get(handles.axes1,'CurrentPoint');     % button up detected
    point1 = point1(1,1:2);                          % extract x and y
    point2 = point2(1,1:2);
end
Xmin = round(min(point1(1,1), point2(1,1)));
Xmax = round(max(point1(1,1), point2(1,1)));
Ymin = round(min(point1(1,2), point2(1,2)));
Ymax = round(max(point1(1,2), point2(1,2)));

Xmin = min(max(Xmin,1),Sx);
Xmax = min(max(Xmax,Xmin+1),Sx);
Ymin = min(max(Ymin,1),Sy);
Ymax = min(max(Ymax,Ymin+1),Sy);

set(handles.edit_Xmin,'String',num2str(Xmin))
set(handles.edit_Xmax,'String',num2str(Xmax))
set(handles.edit_Ymin,'String',num2str(Ymin))
set(handles.edit_Ymax,'String',num2str(Ymax))
draw(hObject, eventdata, handles)

% --- Executes on button press in pushbutton_resetZoom.
function pushbutton_resetZoom_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
set(handles.edit_Xmin,'String',num2str(1))
set(handles.edit_Xmax,'String',num2str(size(Data.I,1)))
set(handles.edit_Ymin,'String',num2str(1))
set(handles.edit_Ymax,'String',num2str(size(Data.I,2)))
draw(hObject, eventdata, handles)


function edit_deltaX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_deltaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_deltaX as text
%        str2double(get(hObject,'String')) returns contents of edit_deltaX as a double
draw(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_deltaX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_deltaX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_deltaY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_deltaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_deltaY as text
%        str2double(get(hObject,'String')) returns contents of edit_deltaY as a double
draw(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_deltaY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_deltaY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_flasePositives.
function radiobutton_flasePositives_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_flasePositives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_flasePositives

if get(handles.radiobutton_flasePositives,'Value')
    set(handles.radiobutton_falseNegatives,'Value',0)
    set(handles.radiobutton_Questionable,'Value',0)
end


% --- Executes on button press in radiobutton_falseNegatives.
function radiobutton_falseNegatives_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_falseNegatives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_falseNegatives

if get(handles.radiobutton_falseNegatives,'Value')
    set(handles.radiobutton_flasePositives,'Value',0)
    set(handles.radiobutton_Questionable,'Value',0)
end


% --- Executes on button press in radiobutton_Questionable.
function radiobutton_Questionable_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Questionable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton_Questionable,'Value')
    set(handles.radiobutton_falseNegatives,'Value',0)
    set(handles.radiobutton_flasePositives,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_Questionable


% --------------------------------------------------------------------
function menu_validateStalls_Callback(hObject, eventdata, handles)
% hObject    handle to menu_validateStalls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
if strcmpi(get(handles.menu_validateStalls,'Checked'), 'off')
    if isfield(Data,'GTStallingMatrix') && isfield(Data,'AutoStallingMatrix')
        set(handles.menu_validateStalls,'Checked','on')
        set(handles.uipanel_validationPanel,'Visible','on')
        set(handles.uipanel_GroundTruthControl,'Visible','on')
    end
    if isfield(Data,'seg')
        if isfield(Data.seg,'LRimage')
            updateStallandFrame(handles, 'NextorCurrent');
            draw(hObject, eventdata, handles);
            makeSegLengthHistogram(handles);
        end
    end
else
    set(handles.menu_validateStalls,'Checked','off')
    set(handles.uipanel_validationPanel,'Visible','off')
    set(handles.uipanel_GroundTruthControl,'Visible','off')
    delete(handles.axes4.Children)
end



function makeSegLengthHistogram(handles)
    global Data
    for i = 1:length(Data.seg)
        pixelLength(1,i) = size(Data.seg(i).LRimage,1);
    end
    delete(handles.axes4.Children)
    hold on
    pixelhist = histogram(handles.axes4,pixelLength,length(pixelLength)-1);    
    handles.axes4.XAxis.Label.String = 'Pixel Length';
    handles.axes4.XAxis.FontWeight = 'Bold';
    handles.axes4.XAxis.TickValues = 0:5:max(pixelLength)+2;
    handles.axes4.YAxis.Label.String = 'Num of Capillary';
    handles.axes4.YAxis.FontWeight = 'Bold';
    handles.text_NumofFilteredValue.String = num2str(sum(pixelLength <= str2num(handles.edit_pixelCutofflengthValue.String)));
    grid on
    Data.segAnalysis.cutoff = str2num(handles.edit_pixelCutofflengthValue.String);
    Data.segAnalysis.capNum = ind2sub(size(pixelLength),find(pixelLength <=Data.segAnalysis.cutoff));
    
    for j = 1:length(Data.seg)
        for i = 1:350
            oneRow = Data.seg(j).LRimage(:,i);
            Data.segAnalysis.AveIntensity(j,i) = mean(oneRow);
            Data.segAnalysis.AveCOV(j,i) = mean(oneRow)/std(oneRow);
        end
        Data.segAnalysis.DatasetAveInt(j,1) = mean(Data.segAnalysis.AveIntensity(j,:)/max(Data.segAnalysis.AveIntensity(j,:)));
        Data.segAnalysis.DatasetAveCOV(j,1) = mean(Data.segAnalysis.AveCOV(j,:)/max(Data.segAnalysis.AveCOV(j,:)));
    end
    

function updateStallandFrame(handles, mov_dir)

global Data
seg_no = str2double(get(handles.edit_segno,'string'));
frame_no = str2double(get(handles.edit_volnumber,'string'));

if get(handles.radiobutton_flasePositives,'Value')
    possible_idx = find(Data.StallingMatrix' == 0 & Data.AutoStallingMatrix' == 1 & Data.ValidationFlag' == 0);
elseif get(handles.radiobutton_falseNegatives,'Value')
    possible_idx = find(Data.StallingMatrix' == 1 & Data.AutoStallingMatrix' == 0 & Data.ValidationFlag' == 0);
elseif get(handles.radiobutton_Questionable,'Value')
    possible_idx = find(Data.GTStallingMatrix' == 2 & Data.ValidationFlag' == 1);
end
if ~isempty(possible_idx)
    idx_to_move = [];
    current_pos_sub = sub2ind(size(Data.StallingMatrix'),frame_no,seg_no);
    if strcmp(mov_dir,'Next')
        idx = find(possible_idx > current_pos_sub);
        if ~isempty(idx)
            idx_to_move = possible_idx(idx(1));
        end
    elseif strcmp(mov_dir,'NextorCurrent')
        idx = find(possible_idx >= current_pos_sub);
        if ~isempty(idx)
            idx_to_move = possible_idx(idx(1));
        end
    elseif strcmp(mov_dir,'Prev')
        idx = find(possible_idx < current_pos_sub);
        if ~isempty(idx)
            idx_to_move = possible_idx(idx(end));
        end
    end
    if ~isempty(idx_to_move)
%         idx_to_move = possible_idx(idx(1));
        seg_no_before = seg_no;
        [frame_no, seg_no] = ind2sub(size(Data.StallingMatrix'),idx_to_move);
        set(handles.edit_segno,'string',num2str(seg_no))
        set(handles.edit_volnumber,'string',num2str(frame_no))
        set(handles.slider_movedata,'Value',frame_no);
        (handles);
        draw([], [], handles);
    end
end


% --- Executes on button press in pushbutton_prevStallforVerification.
function pushbutton_prevStallforVerification_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_prevStallforVerification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateStallandFrame(handles, 'Prev')


% --- Executes on button press in pushbutton_nextStallforVerification.
function pushbutton_nextStallforVerification_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nextStallforVerification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateStallandFrame(handles, 'Next')


% --- Executes on button press in pushbutton_isaStall.
function pushbutton_isaStall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_isaStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
if strcmpi(get(handles.menu_validateStalls,'Checked'), 'on')
    seg_no = str2double(get(handles.edit_segno,'string'));
    frame_no = str2double(get(handles.edit_volnumber,'string'));
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,frame_no) = 1;
    end
    if isfield(Data,'ValidationFlag')
        Data.ValidationFlag(seg_no,frame_no) = 1;
    end
    draw(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton_notaStall.
function pushbutton_notaStall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_notaStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
if strcmpi(get(handles.menu_validateStalls,'Checked'), 'on')
    seg_no = str2double(get(handles.edit_segno,'string'));
    frame_no = str2double(get(handles.edit_volnumber,'string'));
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,frame_no) = 0;
    end
    if isfield(Data,'ValidationFlag')
        Data.ValidationFlag(seg_no,frame_no) = 1;
    end
    draw(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton_Questionable.
function pushbutton_Questionable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Questionable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
if strcmpi(get(handles.menu_validateStalls,'Checked'), 'on')
    seg_no = str2double(get(handles.edit_segno,'string'));
    frame_no = str2double(get(handles.edit_volnumber,'string'));
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,frame_no) = 2;
    end
    if isfield(Data,'ValidationFlag')
        Data.ValidationFlag(seg_no,frame_no) = 1;
    end
    draw(hObject, eventdata, handles)
end


% --- Executes on button press in checkbox_Axe1displayCurrentSegment.
function checkbox_Axe1displayCurrentSegment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Axe1displayCurrentSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Axe1displayCurrentSegment
draw(hObject, eventdata, handles)



function edit_pixelCutofflengthValue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pixelCutofflengthValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pixelCutofflengthValue as text
%        str2double(get(hObject,'String')) returns contents of edit_pixelCutofflengthValue as a double
makeSegLengthHistogram(handles)


% --- Executes during object creation, after setting all properties.
function edit_pixelCutofflengthValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pixelCutofflengthValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_SkipShortSeg.
function pushbutton_SkipSeg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SkipShortSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
if isfield(Data,'GTStallingMatrix')
    Data.GTStallingMatrix(seg_no,:) = 3;
    Data.ValidationFlag(seg_no,:) = 1;
end
disp("Capillary " + num2str(seg_no) + " has been skipped")
updateStallandFrame(handles, 'Next')


% --- Executes on button press in pushbutton_unSkipShortSeg.
function pushbutton_unSkipSeg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_unSkipShortSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
if isfield(Data,'GTStallingMatrix')
    Data.GTStallingMatrix(seg_no,:) = 0;
    Data.ValidationFlag(seg_no,:) = 0;
    
    % When unskip, check for human and autostall matrix again, write
    % GTstalling and Validationflag to be 1 if match
    Manual_StallIndex = Data.StallingMatrix(seg_no,:) == 1;
    Auto_StallIndex = Data.AutoStallingMatrix(seg_no,:) == 1;
    
    % when munual matches auto = 1, GT = 1, flag = 1
    % when manula matches auto = 0, GT = 0, flag = 1
    Data.ValidationFlag(seg_no,:) = Manual_StallIndex == Auto_StallIndex;
    Data.GTStallingMatrix(seg_no,:) = ...
        Data.StallingMatrix(seg_no,:) == 1 &...
        Data.AutoStallingMatrix(seg_no,:) == 1;
    
   disp("Capillary " + num2str(seg_no) + " has been restored to original!")
end
updateStallandFrame(handles, 'Prev')

% --- Executes during object creation, after setting all properties.
function text_NumofFilteredValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NumofFilteredValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global Data
if isfield(Data,'I')
    KeyPressed = upper(eventdata.Key);
    if isfield(Data,'seg')
        if strcmpi(KeyPressed,handles.edit_PreviousCapKey.String) % Previous Cap
            pushbutton_prevseg_Callback(hObject, eventdata, handles);
            disp('Previous Capillary')
        end

        if strcmpi(KeyPressed,handles.edit_NextCapKey.String) % Next Cap
            pushbutton_nextseg_Callback(hObject, eventdata, handles);
            disp('Next Capillary')
        end
    end
    
    if strcmpi(KeyPressed,handles.edit_PreviousImageKey.String) % last image
        pushbutton_moveleft_Callback(hObject, eventdata, handles);
        disp('Previous Image')
    end

    if strcmpi(KeyPressed,handles.edit_NextImageKey.String) % next image
        pushbutton_moveright_Callback(hObject, eventdata, handles);
        disp('Next Image')
    end
    
    if strcmpi(KeyPressed,handles.edit_IncreaseMaxKey.String) % increase upper throshold of contrast
        handles.edit_MaxI.String = num2str(str2double(handles.edit_MaxI.String) + 100);
        edit_MaxI_Callback(hObject, eventdata, handles);
        disp('Increase Max Threshold')
    end
    
    if strcmpi(KeyPressed,handles.edit_DecreaseMaxKey.String) % decrease upper throshold of contrast
        handles.edit_MaxI.String = num2str(str2double(handles.edit_MaxI.String) - 100);
        edit_MaxI_Callback(hObject, eventdata, handles);
        disp('Decrease Max Threshold')
    end

    % Validation Use
    if ~strcmpi(get(handles.menu_validateStalls,'Checked'), 'off')
        if strcmpi(KeyPressed,handles.edit_PreviousCaseKey.String) % Prev
            pushbutton_prevStallforVerification_Callback(hObject, eventdata, handles);
            disp('Previous Case for validation')
        end

        if strcmpi(KeyPressed,handles.edit_NextCaseKey.String) % Next
            pushbutton_nextStallforVerification_Callback(hObject, eventdata, handles);
            disp('Next Case for validation')
        end
        
        if strcmpi(KeyPressed,handles.edit_MarkStallKey.String) % Stall
            pushbutton_isaStall_Callback(hObject, eventdata, handles)
            disp('Select Stall')
        end
        
        if strcmpi(KeyPressed,handles.edit_MarkNotStallKey.String) % Non-Stall
            pushbutton_notaStall_Callback(hObject, eventdata, handles)
            disp('Select Non-Stall')
        end
        
        if strcmpi(KeyPressed,handles.edit_MarkQuestionableKey.String) % Questionable
            pushbutton_Questionable_Callback(hObject, eventdata, handles)
            disp('Select Questionable')
        end
    end
else
    disp('void')
    return
end


% --- Executes on button press in checkbox_IntensityDisp.
function checkbox_IntensityDisp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_IntensityDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_IntensityDisp
global Data
handles.checkbox_COVDisp.Value = 0;
handles.checkbox_IntensityDisp.Value = 0; 

% --- Executes on button press in checkbox_COVDisp.
function checkbox_COVDisp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_COVDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_COVDisp
global Data
handles.checkbox_COVDisp.Value = 0;
handles.checkbox_IntensityDisp.Value = 0; 


% --- Executes on button press in checkbox_Axe2displayCurrentSegment.
function checkbox_Axe2displayCurrentSegment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Axe2displayCurrentSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Axe2displayCurrentSegment
draw(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_Restore.
function pushbutton_Restore_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
start_frame = str2double(handles.edit_GTStartFrame.String);
end_frame = str2double(handles.edit_GTEndFrame.String);
if end_frame <= start_frame
    error('Ending frame is before Starting frame! Please enter the correct range!')
else
    if isfield(Data,'ValidationFlag')
        Data.GTStallingMatrix(seg_no,start_frame:end_frame) = 0;
        Data.ValidationFlag(seg_no,start_frame:end_frame) = 0;
        
        % When unskip, check for human and autostall matrix again, write
        % GTstalling and Validationflag to be 1 if match
        Manual_StallIndex = Data.StallingMatrix(seg_no,start_frame:end_frame) == 1;
        Auto_StallIndex = Data.AutoStallingMatrix(seg_no,start_frame:end_frame) == 1;
        
        % when munual matches auto = 1, GT = 1, flag = 1
        % when manula matches auto = 0, GT = 0, flag = 1
        Data.ValidationFlag(seg_no,start_frame:end_frame) = Manual_StallIndex == Auto_StallIndex;
        Data.GTStallingMatrix(seg_no,start_frame:end_frame) = ...
            Data.StallingMatrix(seg_no,start_frame:end_frame) == 1 &...
            Data.AutoStallingMatrix(seg_no,start_frame:end_frame) == 1;
        
        disp("Frame " + num2str(start_frame) + " to Frame " + num2str(end_frame) ...
            + " on Capillary " + num2str(seg_no) + " has been restored to original!");
    end
    updateStallandFrame(handles, 'Prev')
end


% --- Executes on button press in pushbutton_setNotStall.
function pushbutton_setNotStall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setNotStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
start_frame = str2double(handles.edit_GTStartFrame.String);
end_frame = str2double(handles.edit_GTEndFrame.String);
if end_frame <= start_frame
    error('Ending frame is before Starting frame! Please enter the correct range!')
else
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,start_frame:end_frame) = 0;
        Data.ValidationFlag(seg_no,start_frame:end_frame) = 1;
        disp("Frame " + num2str(start_frame) + " to Frame " + num2str(end_frame) +...
            " on Capillary " + num2str(seg_no) + " do not stall!");
    end
    updateStallandFrame(handles, 'Next')
end



function edit_GTStartFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GTStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_GTStartFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_GTStartFrame as a double


% --- Executes during object creation, after setting all properties.
function edit_GTStartFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GTStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_GTEndFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_GTEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_GTEndFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_GTEndFrame as a double


% --- Executes during object creation, after setting all properties.
function edit_GTEndFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GTEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_setGTStall.
function pushbutton_setGTStall_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setGTStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
start_frame = str2double(handles.edit_GTStartFrame.String);
end_frame = str2double(handles.edit_GTEndFrame.String);
if end_frame <= start_frame
    error('Ending frame is before Starting frame! Please enter the correct range!')
else
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,start_frame:end_frame) = 1;
        Data.ValidationFlag(seg_no,start_frame:end_frame) = 1;
        disp("Frame " + num2str(start_frame) + " to Frame " + num2str(end_frame) +...
            " on Capillary " + num2str(seg_no) + " stalls!");
    end
    updateStallandFrame(handles, 'Next')
end


% --- Executes during object creation, after setting all properties.
function edit_VFlag_Start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_VFlag_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_VFlag_End_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_VFlag_End (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function figure1_WindowKeyPressFcn(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_setGTQuestionable.
function pushbutton_setGTQuestionable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setGTQuestionable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
start_frame = str2double(handles.edit_GTStartFrame.String);
end_frame = str2double(handles.edit_GTEndFrame.String);
if end_frame <= start_frame
    error('Ending frame is before Starting frame! Please enter the correct range!')
else
    if isfield(Data,'GTStallingMatrix')
        Data.GTStallingMatrix(seg_no,start_frame:end_frame) = 2;
        Data.ValidationFlag(seg_no,start_frame:end_frame) = 1;
        disp("Frame " + num2str(start_frame) + " to Frame " + num2str(end_frame) +...
            " on Capillary " + num2str(seg_no) + " is questionable!");
    end
    updateStallandFrame(handles, 'Next')
end



function edit_PreviousCap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousCap as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousCap as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousCap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextCap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextCap as text
%        str2double(get(hObject,'String')) returns contents of edit_NextCap as a double


% --- Executes during object creation, after setting all properties.
function edit_NextCap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousCap as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousCap as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextCap as text
%        str2double(get(hObject,'String')) returns contents of edit_NextCap as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextCap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PreviousImage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousImage as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousImage as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextImage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextImage as text
%        str2double(get(hObject,'String')) returns contents of edit_NextImage as a double


% --- Executes during object creation, after setting all properties.
function edit_NextImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IncreaseMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IncreaseMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IncreaseMax as text
%        str2double(get(hObject,'String')) returns contents of edit_IncreaseMax as a double


% --- Executes during object creation, after setting all properties.
function edit_IncreaseMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IncreaseMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DecreaseMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DecreaseMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DecreaseMax as text
%        str2double(get(hObject,'String')) returns contents of edit_DecreaseMax as a double


% --- Executes during object creation, after setting all properties.
function edit_DecreaseMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DecreaseMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PreviousCase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousCase as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousCase as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousCase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextCase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextCase as text
%        str2double(get(hObject,'String')) returns contents of edit_NextCase as a double


% --- Executes during object creation, after setting all properties.
function edit_NextCase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkStall_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkStall as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkStall as a double


% --- Executes during object creation, after setting all properties.
function edit_MarkStall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkNotStall_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkNotStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkNotStall as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkNotStall as a double


% --- Executes during object creation, after setting all properties.
function edit_MarkNotStall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkNotStall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkQuestionable_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkQuestionable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkQuestionable as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkQuestionable as a double



% --- Executes during object creation, after setting all properties.
function edit_MarkQuestionable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkQuestionable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


    
    
    



function edit_PreviousCapKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCapKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousCapKey as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousCapKey as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousCapKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCapKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextCapKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextCapKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextCapKey as text
%        str2double(get(hObject,'String')) returns contents of edit_NextCapKey as a double


% --- Executes during object creation, after setting all properties.
function edit_NextCapKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextCapKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PreviousImageKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousImageKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousImageKey as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousImageKey as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousImageKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousImageKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextImageKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextImageKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextImageKey as text
%        str2double(get(hObject,'String')) returns contents of edit_NextImageKey as a double


% --- Executes during object creation, after setting all properties.
function edit_NextImageKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextImageKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IncreaseMaxKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IncreaseMaxKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IncreaseMaxKey as text
%        str2double(get(hObject,'String')) returns contents of edit_IncreaseMaxKey as a double


% --- Executes during object creation, after setting all properties.
function edit_IncreaseMaxKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IncreaseMaxKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DecreaseMaxKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DecreaseMaxKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DecreaseMaxKey as text
%        str2double(get(hObject,'String')) returns contents of edit_DecreaseMaxKey as a double


% --- Executes during object creation, after setting all properties.
function edit_DecreaseMaxKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DecreaseMaxKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PreviousCaseKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCaseKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PreviousCaseKey as text
%        str2double(get(hObject,'String')) returns contents of edit_PreviousCaseKey as a double


% --- Executes during object creation, after setting all properties.
function edit_PreviousCaseKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PreviousCaseKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NextCaseKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NextCaseKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NextCaseKey as text
%        str2double(get(hObject,'String')) returns contents of edit_NextCaseKey as a double


% --- Executes during object creation, after setting all properties.
function edit_NextCaseKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NextCaseKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkStallKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkStallKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkStallKey as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkStallKey as a double


% --- Executes during object creation, after setting all properties.
function edit_MarkStallKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkStallKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkNotStallKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkNotStallKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkNotStallKey as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkNotStallKey as a double


% --- Executes during object creation, after setting all properties.
function edit_MarkNotStallKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkNotStallKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MarkQuestionableKey_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MarkQuestionableKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MarkQuestionableKey as text
%        str2double(get(hObject,'String')) returns contents of edit_MarkQuestionableKey as a double


% --- Executes during object creation, after setting all properties.
function edit_MarkQuestionableKey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MarkQuestionableKey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_RestoreKeyDefault.
function pushbutton_RestoreKeyDefault_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RestoreKeyDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.edit_PreviousCapKey.String = 'A';
    handles.edit_NextCapKey.String = 'D';
    handles.edit_PreviousImageKey.String = 'LEFTARROW';
    handles.edit_NextImageKey.String = 'RIGHTARROW';
    handles.edit_IncreaseMaxKey.String = 'W';
    handles.edit_DecreaseMaxKey.String = 'S';
    handles.edit_PreviousCaseKey.String = 'UPARROW';
    handles.edit_NextCaseKey.String = 'DOWNARROW';
    handles.edit_MarkStallKey.String = 'K';
    handles.edit_MarkNotStallKey.String = 'J';
    handles.edit_MarkQuestionableKey.String = 'L';
    disp("Keyboard Shortcuts are restored!")
    
    