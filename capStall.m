function varargout = capStall(varargin)
% CAPSTALL MATLAB code for capStall.fig
%      CAPSTALL, by itself, creates a new CAPSTALL or raises the existing
%      singleton*.
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
%      *See GUI Options on GUIDE's Tools menu_loadresults.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help capStall

% Last Modified by GUIDE v2.5 29-Oct-2021 12:07:07

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
function menu_loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to menu_loaddata (see GCBO)
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

[filename,pathname] = uigetfile({'*.mat;*.tiff;*.tif'},'Please select the Angiogram Data');
if filename == 0
    return
end
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

% pathname = uigetdir('Slect directory for loading data');
% Data.pathname = pathname;
% files = dir([pathname  '/*_angio.mat']);
% load([pathname '/' files(1).name]);
% [x,y] = size(angio);
% z = length(files);
% I = zeros([x y z]);
% for  u = 1:z
%     load([pathname '/' files(u).name]);
%     I(:,:,u) = angio;
% end
% Data.I = I;
% files = dir([pathname  '/*_angioVolume.mat']);
% load([pathname '/' files(1).name]);
% Data.Volume = angio;
% Vz = size(Data.Volume,1);

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

%set(handles.edit_MIPstartidx,'String',num2str(1));
%set(handles.edit_MIPnofframes,'String',num2str(Vz));

% set(handles.slider_moveinZ,'max',Vz);
% set(handles.slider_moveinZ,'min',1);
% set(handles.slider_moveinZ,'Value',Vz);
% set(handles.slider_moveinZ,'SliderStep',[1/(Vz-1), 10/(Vz-1)]);

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
% startidx = str2double(get(handles.edit_MIPstartidx,'String'));
% endidx = str2double(get(handles.edit_MIPnofframes,'String'));
% set(handles.text2,'String',['to ' num2str(startidx+endidx-1) ' spanning']);
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
set(gcf, 'WindowKeyPressFcn', {@figure_WindowKeyPressFcn, handles},'Interruptible','off','BusyAction','cancel');

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
            else
                msg_display = 'Not a Stall';
                color = 'k';
            end
            set(handles.text_GTstallinfo,'String',msg_display, 'fontsize', 10)
            set(handles.text_GTstallinfo,'ForegroundColor',color)
        end
    end
end

if handles.displayCurrentSegment.Value == 1
    if isfield(Data.seg,'pos')
        hold on
        if isfield(Data.seg(seg_no),'frame_seg_pos')
            plot(Data.seg(seg_no).frame_seg_pos(:,2,frame_no),Data.seg(seg_no).frame_seg_pos(:,1,frame_no),'r.','markersize',16);
        else
            plot(Data.seg(seg_no).pos(:,1),Data.seg(seg_no).pos(:,2),'r.','markersize',16);
        end
        hold off
    else
        handles.displayCurrentSegment.Value = 0;
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
        if handles.displayCurrentSegment.Value == 1
            if isfield(Data.seg,'pos')
                hold on
                if isfield(Data.seg(seg_no),'frame_seg_pos')
                    plot(Data.seg(seg_no).frame_seg_pos(:,2,frame_no),Data.seg(seg_no).frame_seg_pos(:,1,frame_no),'r.','markersize',16);
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


% h2 = imagesc(squeeze(max(Data.Volume(startidx:startidx+endidx-1,:,:),[],1)),[MinI MaxI]);
% hold on
% if (get(handles.radiobutton_showseg,'Value') == 1)
%     if isfield(Data,'seg')
%         jj = str2double(get(handles.edit_segno,'string'));
%         plot(Data.seg(jj).pos(:,1),Data.seg(jj).pos(:,2),'r.','markersize',16);
%     end
% else
%     if isfield(Data,'Cap')
%         for u = 1:size(Data.Cap,1)
%             hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','g','FontSize',10);
%             set(hpt,'ButtonDownFcn', sprintf('Cap_Stall_deletept(%d)',u) );
%         end
%     end
%     if isfield(Data,'pts')
%         for u = 1:size(Data.pts,1)
%             hpt1 = text(Data.pts(u,2),Data.pts(u,1),num2str(u),'Color','m','FontSize',10);
%             set(hpt1,'ButtonDownFcn', sprintf('Cap_Stall_deleteZpt(%d)',u) );
%         end
%     end
% end    
%     if isfield(Data,'pts1')
%         for u = 1:size(Data.pts1,1)
%             hpt1 = text(Data.pts1(u,2),Data.pts1(u,1),num2str(u),'Color','m','FontSize',10);
%             set(hpt1,'ButtonDownFcn', sprintf('Cap_Stall_deleteZpt(%d)',3*u-2) );
%         end
%     end
%     if isfield(Data,'pts2')
%         for u = 1:size(Data.pts2,1)
%             hpt2 = text(Data.pts2(u,2),Data.pts2(u,1),num2str(u),'Color','m','FontSize',10);
%             set(hpt2,'ButtonDownFcn', sprintf('Cap_Stall_deleteZpt(%d)',3*u-1) );
%         end
%     end
%     if isfield(Data,'pts3')
%         for u = 1:size(Data.pts3,1)
%             hpt3 = text(Data.pts3(u,2),Data.pts3(u,1),num2str(u),'Color','b','FontSize',10);
%             set(hpt3,'ButtonDownFcn', sprintf('Cap_Stall_deleteZpt(%d)',3*u) );
%         end
%     end

% hold off
% axis image;
% set(h2, 'ButtonDownFcn', {@axes2_ButtonDown, handles});

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
    line(x,y, 'LineWidth', 1, 'Color','white');
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
    end
    hold off
    handles.SegLengthDisplay.String = [ num2str(x(2)-x(1)) ' Pixels' ];
    set(handles.axes3, 'ButtonDownFcn', {@axes3_ButtonDownFcn, handles});
    
    
%     jj = str2double(get(handles.edit_segno,'string'));
%     axes(handles.axes3)
%     hp = plot(Data.Int_ts(jj,:));
%     hold on
%     y = ylim;
%     x = [ii ii];
%     line(x,y);
%     if isfield(Data,'StallingPts')
%         for uu = 1:length(Data.StallingPts(jj).pos)
%             xpt = Data.StallingPts(jj).pos(uu);
%             ypt = Data.Int_ts(jj,xpt);
%             text(xpt,ypt,'*');
%         end
%     end
%     if isfield(Data,'StallingMatrix')
%         xidx = find(Data.StallingMatrix(jj,:) == 1);
%         yidx = Data.Int_ts(jj,xidx);
%         text(xidx,yidx,'*');
%     end
%     hold off
%     ylim(y);
%     set(handles.axes3, 'ButtonDownFcn', {@axes3_ButtonDownFcn, handles});
end


Data.handles =  handles;
Data.hObject =  hObject;
Data.eventdata = eventdata;


%  if isfield(Data,'sliderobject')
%     uicontrol(Data.sliderobject);
% end
% set(handles,'Keypressfcn',@keypressed);
function figure_WindowKeyPressFcn(hObject, eventdata, handles)

global Data
% 
% if ~strcmpi(get(gco,'Tag'),'figure1')
%     if strcmpi(get(gco,'style'),'edit')
%         return;
%     end
% end
% 
% keyPressed = eventdata.Key;
% if strcmpi(keyPressed,'space')
%     if isfield(Data,'Int_ts')
%         jj = str2double(get(handles.edit_segno,'string'));
%         ii = str2double(get(handles.edit_volnumber,'string'));
%         if isfield(Data,'StallingMatrix')
%             if(Data.StallingMatrix(jj,ii) == 0)
%                 Data.StallingMatrix(jj,ii) = 1;
%             else
%                 Data.StallingMatrix(jj,ii) = 0;
%             end
%         else
%             StallingMatrix = zeros(length(Data.seg),size(Data.I,3));
%             StallingMatrix(jj,ii) = 1;
%             Data.StallingMatrix = StallingMatrix; 
%         end
%     end
%     draw(hObject, eventdata, handles);
% end
% if isfield(Data,'sliderobject')
%     uicontrol(Data.sliderobject);
% end
 


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
% elseif strcmp(mouseclick,'alt')
%      ii = round(eventdata.IntersectionPoint(1));
%      jj = str2double(get(handles.edit_segno,'string'));
%      if isfield(Data,'StallingPts')
%          if isfield(Data.StallingPts,'pos')
%              Data.StallingPts(jj).pos = [Data.StallingPts(jj).pos; ii];
%          else
%              Data.StallingPts(jj).pos = ii;
%          end
%      else 
%          stallcount = length(Data.seg);
%          Data.StallingPts(stallcount).pos = [];
%          Data.StallingPts(jj).pos = ii;
%      end
% end
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
%     return
% end
% if (axis1pos(1) >=1 && axis1pos(1) <= I_x && axis1pos(2) >=1 && axis1pos(2) <= I_y)
%     MinI = str2double(get(handles.edit_MinI,'String'));
%     MaxI = str2double(get(handles.edit_MaxI,'String'));
%     if eventdata.VerticalScrollCount > 0
%         ii = str2num(get(handles.edit_volnumber,'string'));
%         ii = min(max(ii+1,1),size(Data.I,3));
%         set(handles.edit_volnumber,'string',num2str(ii));
%         set(handles.slider_movedata,'Value',ii);
% %         draw(hObject, eventdata, handles);
%         axes(handles.axes1)
%         colormap('gray');
%         h = imagesc(Data.I(:,:,ii),[MinI MaxI]);
%         hold on
%         if isfield(Data,'Cap')
%             for u = 1:size(Data.Cap,1)
%                 if u == cseg
%                     hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','b','FontSize',10);
%                 else
%                     hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','g','FontSize',10);
%                 end
%                 set(hpt,'ButtonDownFcn', sprintf('Cap_Stall_deletept(%d)',u) );
%             end
%         end
%         hold off
%         axis image;
%         set(h, 'ButtonDownFcn', {@axes1_ButtonDown, handles});
%         set(gcf, 'WindowScrollWheelFcn', {@axes_WindowScrollWheelFcn, handles},'Interruptible','off');
%         set(gcf, 'WindowKeyPressFcn', {@figure_WindowKeyPressFcn, handles},'Interruptible','off');
%     end  
%     if eventdata.VerticalScrollCount < 0
%         ii = str2num(get(handles.edit_volnumber,'string'));
%         ii = min(max(ii-1,1),size(Data.I,3));
%         set(handles.edit_volnumber,'string',num2str(ii));
%         set(handles.slider_movedata,'Value',ii);
% %         draw(hObject, eventdata, handles); 
%         axes(handles.axes1)
%         colormap('gray');
%         h = imagesc(Data.I(:,:,ii),[MinI MaxI]);
%         hold on
%         if isfield(Data,'Cap')
%             for u = 1:size(Data.Cap,1)
%                 if u == cseg
%                     hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','b','FontSize',10);
%                 else
%                     hpt = text(Data.Cap(u,2),Data.Cap(u,1),num2str(u),'Color','g','FontSize',10);
%                 end
%                 set(hpt,'ButtonDownFcn', sprintf('Cap_Stall_deletept(%d)',u) );
%             end
%         end
%         hold off
%         axis image;
%         set(h, 'ButtonDownFcn', {@axes1_ButtonDown, handles});
%         set(gcf, 'WindowScrollWheelFcn', {@axes_WindowScrollWheelFcn, handles},'Interruptible','off');
%         set(gcf, 'WindowKeyPressFcn', {@figure_WindowKeyPressFcn, handles},'Interruptible','off');
%     end  
%     Xmin = str2double(get(handles.edit_Xmin,'String'));
%     Ymin = str2double(get(handles.edit_Ymin,'String'));
%     Xmax = str2double(get(handles.edit_Xmax,'String'));
%     Ymax = str2double(get(handles.edit_Ymax,'String'));
%     xlim([Xmin Xmax])
%     ylim([Ymin Ymax])
%     
%     if isfield(Data,'sliderobject')
%         uicontrol(Data.sliderobject);
%     end
% elseif axis2pos(1) >=1 && axis2pos(1) <= V_x && axis2pos(2) >=1 && axis2pos(2) <= V_y
%     if eventdata.VerticalScrollCount > 0
%         Vz = size(Data.Volume,1);
%         ii = str2double(get(handles.edit_MIPstartidx,'String'));
%         jj = str2double(get(handles.edit_MIPnofframes,'String'));
%         ii = min(max(ii+1,1),Vz);
%         jj = min(max(jj,1),Vz-ii+1);
%         set(handles.edit_MIPstartidx,'String',num2str(ii));
%         set(handles.edit_MIPnofframes,'String',num2str(jj));
%         set(handles.slider_moveinZ,'Value',Vz-ii+1);
%         draw(hObject, eventdata, handles);
%     end
%     if eventdata.VerticalScrollCount < 0
%         Vz = size(Data.Volume,1);
%         ii = str2double(get(handles.edit_MIPstartidx,'String'));
%         jj = str2double(get(handles.edit_MIPnofframes,'String'));
%         ii = min(max(ii-1,1),Vz);
%         jj = min(max(jj,1),Vz-ii+1);
%         set(handles.edit_MIPstartidx,'String',num2str(ii));
%         set(handles.edit_MIPnofframes,'String',num2str(jj));
%         set(handles.slider_moveinZ,'Value',Vz-ii+1);
%         draw(hObject, eventdata, handles);
%     end
%     if isfield(Data,'sliderobjectZ')
%          uicontrol(Data.sliderobjectZ);
%     end
% end
% 
% 
% % pt1 = get(handles.axes1,'Position')
% 





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

% function axes2_ButtonDown(hObject, eventdata, handles)
% 
% global Data
% 
% parent = (get(hObject, 'Parent'));
% pts = round(get(parent, 'CurrentPoint'));
% y = pts(1,1);
% x = pts(1,2);
% z = str2double(get(handles.edit_MIPstartidx,'String'));
% if isfield(Data,'pts')
%     Data.pts = [Data.pts; x y z];
% else
%     Data.pts = [x y z];
% end
% count = 0;
% if isfield(Data,'count')
%     if rem(Data.count,3) == 0
%         if isfield(Data,'pts1')
%             Data.pts1 = [Data.pts1; x y z];
%         else
%             Data.pts1 = [x y z];
%         end
%     end
%     if rem(Data.count,3) == 1
%         if isfield(Data,'pts2')
%             Data.pts2 = [Data.pts2; x y z];
%         else
%             Data.pts2 = [x y z];
%         end
%     end
%     if rem(Data.count,3) == 2
%         if isfield(Data,'pts3')
%             Data.pts3 = [Data.pts3; x y z];
%         else
%             Data.pts3 = [x y z];
%         end
%     end
%     Data.count = Data.count+1;
% else
%     Data.pts1 = [x y z];
%     Data.count = 1;
% end
    
% if isfield(Data,'deletedpts')
%     if(isempty(Data.deletedpts) == 0)
%         if rem(Data.deletedpts(1),3) == 1
%             count = Data.deletedpts(1);
%             Data.deletedpts(1) = [];
%         elseif rem(Data.deletedpts(1),3) == 2
%             count = Data.deletedpts(1);
%             Data.deletedpts(1) = [];
%         elseif rem(Data.deletedpts(1),3) == 0
%             count = Data.deletedpts(1);
%             Data.deletedpts(1) = [];
%         end
%     else
%         if isfield(Data,'count')
%             Data.count = Data.count+1;
%         else
%             Data.count = 1;
%         end
%     end    
% else
%     if isfield(Data,'count')
%         Data.count = Data.count+1;
%     else
%         Data.count = 1;
%     end
% end
% if count == 0
%     if rem(Data.count,3) == 1
%         if isfield(Data,'pts1')
%             Data.pts1 = [Data.pts1; x y z];
%         else
%             Data.pts1 = [x y z];
%         end
%     end
%     if rem(Data.count,3) == 2
%         if isfield(Data,'pts2')
%             Data.pts2 = [Data.pts2; x y z];
%         else
%             Data.pts2 = [x y z];
%         end
%     end
%     if rem(Data.count,3) == 0
%         if isfield(Data,'pts3')
%             Data.pts3 = [Data.pts3; x y z];
%         else
%             Data.pts3 = [x y z];
%         end
%     end
% else 
%     idx = floor((count-1)/3)+1;
%     if rem(count,3) == 1
%         Data.pts1(idx,:) = [x y z];
%     elseif rem(count,3) == 2
%         Data.pts2(idx,:) = [x y z];
%     elseif rem(count,3) == 0
%         Data.pts3(idx,:) = [x y z];
%     end
% end
draw(hObject, eventdata, handles);
if isfield(Data,'sliderobjectZ')
    uicontrol(Data.sliderobjectZ);
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
%     for uuu = 1:size(mean_seg_pos,1)
%         startx = max(mean_seg_pos(uuu,1)-4,1);
%         endx = min(mean_seg_pos(uuu,1)+4,size(mean_seg_pos,1));
%         starty = max(mean_seg_pos(uuu,2)-4,1);
%         endy = min(mean_seg_pos(uuu,2)+4,size(mean_seg_pos,2));
%         mean_mask(startx:endx,starty:endy) = 1;
%     end
%     mean_cropped_image = mean_cropped_image.*mean_mask;
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
% Z = size(Data.I,3);
% Int_ts = zeros(length(seg),Z);
% for kk = 1:length(seg)
%     temp_ts = zeros(size(seg(kk).pos,1),Z);
%     for ll = 1:size(seg(kk).pos,1)
%         temp_ts(ll,:) = squeeze(Data.I(seg(kk).pos(ll,2),seg(kk).pos(ll,1),:));
%     end
%     Int_ts(kk,:) = mean(temp_ts,1);
% end
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
smoothFactor = 6;
corrThresh = 0.7;
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

% create validation flag
if isfield(Data,'ValidationFlag')
    Int_rows = length(Data.seg);
    Stall_rows = size(Data.ValidationFlag,1);
    stall_cols = size(Data.ValidationFlag,2);
    if Int_rows > Stall_rows
       Data.ValidationFlag = [Data.ValidationFlag; zeros([Int_rows-Stall_rows stall_cols])];
    end
else
    Data.ValidationFlag = zeros(length(seg),n_frames);
    idx = find(Data.StallingMatrix == 1 & Data.AutoStallingMatrix==1);
    Data.ValidationFlag(idx) = 1;
    idx = find(Data.StallingMatrix == 0 & Data.AutoStallingMatrix==0);
    Data.ValidationFlag(idx) = 1;
end
draw(hObject, eventdata, handles);
makeSegLengthHistogram(handles);
makeSegQualityAnalysis(handles);

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


%%
% I = Data.Volume;
% Z = size(Data.I,3);
% beta = 100;
% c = 500;
% [k,l,m] = size(I);
% [nz,nx,ny] = size(I);
% L   = zeros(k,l,m);
% Vs  = zeros(k,l,m);
% alpha = 0.25;
% gamma12 = 0.5;
% gamma23 = 0.5;
% T = zeros(k,l,m);
% sigma = [2.5];
%     for i = 1:length(sigma)
%         
% 
%         [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigma(i));
% 
%         %Normalizing the Hessian Matrix
%         %Dxx = i^2*Dxx; Dyy = i^2*Dyy;  Dzz = i^2*Dzz; Dxy = i^2*Dxy;  Dxz = i^2*Dxz; Dyz = i^2*Dyz;
%         
%        
% 
%         [Lambda1,Lambda2,Lambda3,V1,V2,V3] = eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
%         
%         SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'descend');
%         Lambda1 = reshape(SortL(1,:),size(Lambda1));
%         Lambda2 = reshape(SortL(2,:),size(Lambda2));
%         Lambda3 = reshape(SortL(3,:),size(Lambda3));
%         
%         idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 < 0);
%         T(idx ) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1+Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
%         idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 > 0 & Lambda1 < abs(Lambda2)/alpha);
%         T(idx ) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1-alpha*Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
% 
% %         L1 = (2*Lambda1-Lambda2-Lambda3)./(2*sqrt(Lambda1.^2+Lambda2.^2+Lambda3.^2-Lambda1.*Lambda2-Lambda1.*Lambda3-Lambda2.*Lambda3));
% %         L1 = exp(-alpha*(L1-1).^2); 
% %         L1(abs(Lambda1)> abs(Lambda2)) = 0;
% %         L1(Lambda2>0 | Lambda3>0) = 0;
% %         L1 = -L1.*Lambda2;
% %         L = max(L,L1); 
% % 
% %         Ra = abs(Lambda2./Lambda3);
% %         Rb = abs(Lambda1./(Lambda2.*Lambda3));
% %         s = sqrt(Lambda1.^2+Lambda2.^2+Lambda3.^2);
% %         Vs1 = 1-exp(-Ra.^2/(2*alpha)).*exp(-Rb.^2/(2*beta)).*(1-exp(-s.^2/(2*c)));
% %         Vs1(Lambda2>0 | Lambda3>0) = 0;
% %         Vs1(abs(Lambda1) > abs(Lambda2)) = 0;
% %         Vs = max(Vs,Vs1); 
%              [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(T,sigma);
%         
%         [Lambda1,Lambda2,Lambda3,V1,V2,V3] = eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
%         
%         SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'ascend');
%         Lambda1 = reshape(SortL(1,:),size(Lambda1));
%         Lambda2 = reshape(SortL(2,:),size(Lambda2));
%         Lambda3 = reshape(SortL(3,:),size(Lambda3));
%         
%         E = -sigma^2.*Lambda2;
%         E(E<0) = 0;
%     end
% %     T = L;
% T = E;
% 
% T = (T-min(T(:)))*255/(max(T(:))-min(T(:)));
% 
% segm = zeros(size(I));
% %     if isfield(cna,'Dzz')~= 1
%     [Dzz, Dxx, Dyy, Dzx, Dzy, Dxy] = Hessian3D(T,sigma(1));
%     [Lambda1,Lambda2,Lambda3,V1,V2,V3] = eig3volume(Dzz,Dzx,Dzy,Dxx,Dxy,Dyy);
%     cna.Dzz = Dzz; cna.Dxx = Dxx; cna.Dyy = Dyy; cna.Dzx = Dzx; cna.Dzy = Dzy; cna.Dxy = Dxy;
%     Lambda1 = Lambda1; Lambda2 = Lambda2; Lambda3 = Lambda3;
% %     else
% %        Dzz = cna.Dzz; Dxx = cna.Dxx; Dyy = cna.Dyy; Dzx = cna.Dzx; Dzy = cna.Dzy; Dxy = cna.Dxy;
% %        Lambda1 = cna.Lambda1; Lambda2 = cna.Lambda2; Lambda3 = cna.Lambda3; 
% %     end
%     
%     T1 = (T-min(T(:)))/(max(T(:))-min(T(:)));
%     T_seg = zeros(size(T));
%     idx = find(T1 > 0.075);
%     T_seg(idx) = 1;
%     
%     CC = bwconncomp(T_seg);
%     T_segM = T_seg;
%     T_cpt = T_seg*0;
%     for uuu = 1:length(CC.PixelIdxList)
%          if length(CC.PixelIdxList{uuu}) < 100
%              T_segM(CC.PixelIdxList{uuu}) = 0;
%          end
%     end
%     
%     count  = 0;
%     VSidx = find(T_segM == 1);
%     seg_no = 0;
%     for vv = 1:size(Data.pts,1)
%         vv
%         seg_no =seg_no+1;
% 
%         for direction = 1:2
%             z = Data.pts(vv,3);
%             x = Data.pts(vv,1);
%             y = Data.pts(vv,2);
%             xf = x;
%             yf = y;
%             zf = z;
%             coor = zeros([7 3]);
%             coor(1,:) = [z x y];
%             s = 1;
%             BFB = 0;
%             BI = 0;
% 
%             while 1
%                 if s==1
%     %                 [z x y]
%                     M = [ Dzz(z,x,y)  Dzx(z,x,y) Dzy(z,x,y); Dzx(z,x,y) Dxx(z,x,y) Dxy(z,x,y); Dzy(z,x,y)  Dxy(z,x,y) Dyy(z,x,y)];
%                     [V,D] = eig(M);
%                     D = sum(D,1);
%                     [~,id] = min(abs(D));
%                     idx = find(abs(D) ~= min(abs(D)));     VD = V(:,id);
%                      if s > 1
%                         dp = sum(VD.*VP);
%                         if   dp < 0
%                             VD = -VD;
%                         end
%                      end
%                 end
%                 FS = 2;
%                 xc = repmat(-FS:FS,[2*FS+1 1]);
%                 yc = repmat((-FS:FS)',[1 2*FS+1]);
%                 zc = zeros(2*FS+1,2*FS+1);
%                 VC = [1;0;0];
%                 costheta = (VC'*VD)/(norm(VC)*norm(VD));
%                 sintheta = norm(cross(VC,VD))/(norm(VC)*norm(VD));
%                 C = cross(VC,VD);
%                 C = C/norm(C);
%                 xc = xc(:);
%                 yc = yc(:);
%                 zc = zc(:);
% 
%                 xx = (C(2)^2*xc+C(3)*C(3)*yc+C(2)*C(1)*zc)*(1-costheta)+xc*costheta+(C(3)*zc-C(1)*yc)*sintheta;
%                 yy = (C(2)*C(3)*xc+C(3)^2*yc+C(3)*C(1)*zc)*(1-costheta)+yc*costheta+(C(1)*xc-C(2)*zc)*sintheta;
%                 zz= (C(3)*C(1)*xc+C(3)*C(1)*yc+C(3)^2*zc)*(1-costheta)+zc*costheta+(C(2)*yc-C(3)*xc)*sintheta;
% 
%                 xx = xx+xf;
%                 yy = yy+yf;
%                 zz = zz+zf;
%                 img = zeros(size(xc(:)));
% 
%                 for u = 1:size(xc(:))
%                       img(u) = T(min(max(round(zz(u)),1),nz),min(max(round((xx(u))),1),nx),min(max(round(yy(u)),1),ny));
%                 %     img(u) = (cna.angio(ceil(zz(u)),ceil(xx(u)),ceil(yy(u)))+cna.angio(floor(zz(u)),floor(xx(u)),floor(yy(u)))+...
%                 %              cna.angio(ceil(zz(u)),ceil(xx(u)),floor(yy(u)))+cna.angio(ceil(zz(u)),floor(xx(u)),ceil(yy(u)))+...
%                 %              cna.angio(ceil(zz(u)),floor(xx(u)),floor(yy(u)))+cna.angio(floor(zz(u)),ceil(xx(u)),ceil(yy(u)))+...
%                 %              cna.angio(floor(zz(u)),ceil(xx(u)),floor(yy(u)))+cna.angio(floor(zz(u)),floor(xx(u)),ceil(yy(u))))/8;
%                 end
% 
%                 img = reshape(img,[2*FS+1 2*FS+1]);
% 
%                 [~,C1idx] = max(squeeze(mean(img,2)));
%                 [~,C2idx] = max(squeeze(mean(img,1)));
%                 uidx = sub2ind(size(img),C1idx,C2idx);
% 
%                 Cz(s) = min(max(round(zz(uidx)),1),nz);
%                 Cx(s) = min(max(round(xx(uidx)),1),nx);
%                 Cy(s) = min(max(round(yy(uidx)),1),ny);
%                 if direction == 1
%                     C1z(s) = min(max(round(zz(uidx)),1),nz);
%                     C1x(s) = min(max(round(xx(uidx)),1),nx);
%                     C1y(s) = min(max(round(yy(uidx)),1),ny);
%                 end
% 
%                 segm(Cz(s),Cx(s),Cy(s)) = 1;
%                 ss= 7;
%                 surmat = T(min(max(Cz(s)-ss,1),nz):min(max(Cz(s)+ss,1),nz),min(max(Cx(s)-ss,1),nx):min(max(Cx(s)+ss,1),nx),min(max(Cy(s)-ss,1),ny):min(max(Cy(s)+ss,1),ny)) ;
%     %            Vidx = find(surmat > 0.1*img(C1idx,C2idx)); 
%                 Vidx = find(surmat~=0 | surmat==0);
%                 [Vsub1,Vsub2,Vsub3] = ind2sub(size(surmat),Vidx);
%     %            shiftmat = [Cz(s) Cx(s) Cy(s); Cz(s) Cx(s) Cy(s); Cz(s) Cx(s) Cy(s)];
%                 shiftmat = [repmat(Cz(s),[length(Vsub1) 1]) repmat(Cx(s),[length(Vsub2) 1]) repmat(Cy(s),[length(Vsub3) 1])];
% 
%                 Vsub = [Vsub1 Vsub2 Vsub3]-ss*ones(size([Vsub1 Vsub2 Vsub3]))+shiftmat;
%                tidx = find(Vsub(:,1) < 1 |Vsub(:,1) > nz | Vsub(:,2) < 1 | Vsub(:,2) > nx | Vsub(:,3) < 1 |Vsub(:,3) > ny);
%                Vsub(tidx,:,:) = [];
%                VFidx = sub2ind(size(T),Vsub(:,1),Vsub(:,2),Vsub(:,3));
%     %            Vf1idx = sub2idx(size(T));
%                idx_intersect = intersect(VFidx,VSidx);
%                idxvaluept = sub2ind(size(T),z,x,y);
%     %            
%                T_segM([idxvaluept;idx_intersect]) = seg_no;
%                VSidx = setdiff(VSidx,[idxvaluept;idx_intersect]);
%                if s == 1
%                    if direction == 1
%                        mask = [idxvaluept;idx_intersect];
%                    else 
%                        mask = [mask; idxvaluept;idx_intersect];
%                    end
%                else
%                    mask = [mask; idxvaluept;idx_intersect];
%                end
% 
%                vidx = find(img > 0.3*img(C1idx,C2idx));
%                sumArea(s) = length(vidx);
%                sumT(s) = img(C1idx, C2idx);
% 
%                step = 2;
%         %        noRecenter = 0;
%                if s>1
%                     if sum([Cz(s)-Cz(s-1) Cx(s)-Cx(s-1) Cy(s)-Cy(s-1)])==0
%     %                    disp( sprintf('Recenter did not move. prev_step = %d. ',prev_step) )
%                        step = prev_step + 1;
%         %                noRecenter = 1;
%                    end
% 
%         %            if VP(1) >= 1/(step+1) && Cz(s) <= Cz(s-1);
%         %               
%         %                    step = prev_step+1;
%         %                   
%         %            elseif  VP(1) <= -1/(step+1) && Cz(s) >= Cz(s-1);
%         %                     step = prev_step+1;
%         %                  
%         %            elseif VP(2) >= 1/(step+1) && Cx(s) <= Cx(s-1);
%         %                      step = prev_step+1;
%         %                    
%         %             elseif VP(2) <= -1/(step+1) && Cx(s) >= Cx(s-1);
%         %                     step = prev_step+1;
%         %                    
%         %             elseif VP(3) >= 1/(step+1)&& Cy(s) <= Cy(s-1);
%         %                      step = prev_step+1;
%         %                     
%         %             elseif  VP(3) <= -1/(step+1) && Cy(s) >= Cy(s-1);
%         %                     step = prev_step+1;
%         %                     
%         %             end
%                end
% 
%                if step > 4
%                    BFB = BFB+1;
%                    break;
%                end
%                prev_step = step;
% 
%                  M = [ Dzz(Cz(s),Cx(s),Cy(s))  Dzx(Cz(s),Cx(s),Cy(s)) Dzy(Cz(s),Cx(s),Cy(s)); Dzx(Cz(s),Cx(s),Cy(s)) Dxx(Cz(s),Cx(s),Cy(s)) Dxy(Cz(s),Cx(s),Cy(s)); Dzy(Cz(s),Cx(s),Cy(s))  Dxy(Cz(s),Cx(s),Cy(s)) Dyy(Cz(s),Cx(s),Cy(s))];
%                    [V,D] = eig(M);
%                    D = sum(D,1);
%                    [~,id] = min(abs(D));
%                    idx = find(abs(D) ~= min(abs(D)));     VD = V(:,id);
%                    if s== 1
%                        if direction == 2
%                            VD = -VD;
%                        end
%                    end
%                    if s > 1
%                         dp = sum(VD.*VP);
%                         if  dp < 0
%                             VD = -VD;
%                         end
%                    end
% 
%                    if s>1
%                    VPP = VP;
%                end
%                VP = VD;
%                zf = (Cz(s)+step*VD(1)); z = round(zf);
%                xf = (Cx(s)+step*VD(2)); x = round(xf);
%                yf = (Cy(s)+step*VD(3)); y = round(yf);
%                coor(7,:) = coor(6,:);
%                coor(6,:) = coor(5,:);
%                coor(5,:) = coor(4,:);
%                coor(4,:) = coor(3,:);
%                coor(3,:) = coor(2,:);
%                coor(2,:) = coor(1,:);
%                coor(1,:) = [z x y];
% 
%                if z<1 || z>nz || x<1 || x>nx || y<1 || y>ny
% 
%     %                disp('out of bounds');
%                    break;
%                end
% 
%                if coor(1,:) == coor(2,:) | coor(1,:) == coor(3,:) | coor(1,:) == coor(4,:) | coor(1,:) == coor(5,:) | coor(1,:) == coor(6,:) | coor(1,:) == coor(7,:)
%                    BFB = BFB+1;
%                    break;
%                end
% 
%                lambdasC(:,s) =  [Lambda1(Cz(s),Cx(s),Cy(s)); Lambda2(Cz(s),Cx(s),Cy(s)); Lambda3(Cz(s),Cx(s),Cy(s))];
%                if direction ==1
%                    lambdasC1(:,s) =  [Lambda1(Cz(s),Cx(s),Cy(s)); Lambda2(Cz(s),Cx(s),Cy(s)); Lambda3(Cz(s),Cx(s),Cy(s))];
%                end
%                lambdasN(:,s) =  [Lambda1(z,x,y); Lambda2(z,x,y); Lambda3(z,x,y)];
% 
% %                dist1 = sqrt((z-Data.pts1(vv,3))^2+(x-Data.pts1(vv,1))^2+(y-Data.pts1(vv,2))^2);
% %                dist2 = sqrt((z-Data.pts2(vv,3))^2+(x-Data.pts2(vv,1))^2+(y-Data.pts2(vv,2))^2);
% %                if dist1 <= 5 || dist2 <= 5 || s>40
% %     %                disp(sprintf('sumT(%d) = %.3f is <= 10', s, sumT(s)));
% %                    BI = BI+1;
% %                    break;
% %                end
%                 if s > 10
%                     break;
%                 end
%                s = s+1;
%               
% 
%             end  
%             if direction ==1
%                 seg(seg_no).pos = [C1z(end:-1:1)' C1x(end:-1:1)' C1y(end:-1:1)'];
%             elseif direction ==2
%                 seg(seg_no).pos = [seg(seg_no).pos; Cz(1:s)' Cx(1:s)' Cy(1:s)'];
%             end
%             C1z = []; Cz =[]; C1x = []; Cx =[]; C1y = []; Cy =[];
%         end
%         seg(seg_no).mask = mask;
%     end
%     Data.seg = seg;
%     
%     figure; colormap('gray');
%     imagesc(log(squeeze(max(I,[],1)))); 
%     hold on
%     for uu = 1:length(seg)
%         plot(seg(uu).pos(:,3),seg(uu).pos(:,2),'r.','markersize',16);
%     end
%     hold off
%     
%    Int_ts = zeros(length(seg),Z);
%    slashpos = max(strfind(Data.pathname,filesep));
%    datapath = Data.pathname(1:slashpos-1);
%    files = dir([datapath '/*_angio.mat']);
%    underscorepos = strfind(Data.pathname,'_');
%    underscorepos = underscorepos(end-1:end);
%    startidx = str2num(Data.pathname(underscorepos(1)+1:underscorepos(2)-1));
%    endidx = str2num(Data.pathname(underscorepos(2)+1:end));
%        for kk = 1:length(files)
%            kk
%            load([datapath '/' files(kk).name]);
%            angio = angio(startidx:endidx,:,:);
%            for ll = 1:length(seg)
%                Int_ts(ll,kk) = mean(angio(seg(ll).mask));
%            end
%        end
%    
%    Data.Int_ts = Int_ts;
%    if isfield(Data,'StallingMatrix')
%        Int_rows = size(Data.Int_ts,1);
%        Stall_rows = size(Data.StallingMatrix,1);
%        stall_cols = size(Data.StallingMatrix,2);
%        if Int_rows > Stall_rows
%            Data.StallingMatrix = [Data.StallingMatrix; zeros([Int_rows-Stall_rows stall_cols])];
%        end
%    end
%    draw(hObject, eventdata, handles);
    
    
    



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
SkipAllDisplay(handles);
if isfield(Data,'seg')
    if isfield(Data.seg,'LRimage')
        makeSegQualityAnalysis(handles);
        makeDatasetQualityAnalysis(handles);
    end
end
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
SkipAllDisplay(handles);
if isfield(Data,'seg')
    if isfield(Data.seg,'LRimage')
        makeSegQualityAnalysis(handles);
        makeDatasetQualityAnalysis(handles);
    end
end
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
SkipAllDisplay(handles);
if isfield(Data,'seg')
    if isfield(Data.seg,'LRimage')
        makeSegQualityAnalysis(handles);
        makeDatasetQualityAnalysis(handles);
    end
end
draw(hObject, eventdata, handles);


% function update_zoom_with_segno(seg_no, handles)
% 
% global Data
% 
% [Sx, Sy, ~] = size(Data.I);
% seg_pos = Data.seg(seg_no).pos;
% deltaX = str2double(get(handles.edit_deltaX,'String'));
% deltaY = str2double(get(handles.edit_deltaY,'String'));
% Xmin = min(seg_pos(:,1))-deltaX;
% Xmax = max(seg_pos(:,1))+deltaX;
% Ymin = min(seg_pos(:,2))-deltaY;
% Ymax = max(seg_pos(:,2))+deltaY;
% 
% Xmin = min(max(Xmin,1),Sx);
% Xmax = min(max(Xmax,Xmin+1),Sx);
% Ymin = min(max(Ymin,1),Sy);
% Ymax = min(max(Ymax,Ymin+1),Sy);

% set(handles.edit_Xmin,'String',num2str(Xmin))
% set(handles.edit_Xmax,'String',num2str(Xmax))
% set(handles.edit_Ymin,'String',num2str(Ymin))
% set(handles.edit_Ymax,'String',num2str(Ymax))



% --- Executes on button press in radiobutton_showseg.
function radiobutton_showseg_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_showseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_showseg

draw(hObject, eventdata, handles);

% 
% % --- Executes on button press in pushbutton_SaveResults.
% function pushbutton_SaveResults_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_SaveResults (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)





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
function menu_loadResults_Callback(hObject, eventdata, handles)
% hObject    handle to menu_loadResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data;
handles.displayCurrentSegment.Value = 0;
delete(handles.axes3.Children)
delete(handles.axes4.Children)
delete(handles.axes5.Children)
delete(handles.axes6.Children)

[filename, pathname] = uigetfile;
if filename == 0
    return
end
temp_struct = load([pathname filename]);
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
   handles.displayCurrentSegment.Value = 1;
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
    end
    if isfield(Data,'seg')
        if isfield(Data.seg,'LRimage')
            updateStallandFrame(handles, 'NextorCurrent');
            draw(hObject, eventdata, handles);
            makeSegLengthHistogram(handles);
            makeSegQualityAnalysis(handles);
            makeDatasetQualityAnalysis(handles);
        end
    end
%     if get(handles.radiobutton_flasePositives,'Value')
%         seg_no = str2double(get(handles.edit_segno,'string'));
%         frame_no = str2double(get(handles.edit_volnumber,'string'));
%         false_positives_idx = find(Data.StallingMatrix' == 0 & Data.AutoStallingMatrix' == 1 & Data.ValidationFlag' == 0);
%         current_pos_sub = sub2ind(size(Data.StallingMatrix'),frame_no,seg_no);
%         idx = find(false_positives_idx >= current_pos_sub);
%         if ~isempty(idx)
%             idx_to_move = false_positives_idx(idx(1));
%             [frame_no, seg_no] = ind2sub(size(Data.StallingMatrix'),idx_to_move);
%             set(handles.edit_segno,'string',num2str(seg_no))
%             set(handles.edit_volnumber,'string',num2str(frame_no))
%             draw(hObject, eventdata, handles);
%         end
%     elseif get(handles.radiobutton_falseNegatives,'Value')
%         seg_no = str2double(get(handles.edit_segno,'string'));
%         frame_no = str2double(get(handles.edit_volnumber,'string'));
%         false_positives_idx = find(Data.StallingMatrix' == 1 & Data.AutoStallingMatrix' == 0 & Data.ValidationFlag' == 0);
%         current_pos_sub = sub2ind(size(Data.StallingMatrix'),frame_no,seg_no);
%         idx = find(false_positives_idx >= current_pos_sub);
%         if ~isempty(idx)
%             idx_to_move = false_positives_idx(idx(1));
%             [frame_no, seg_no] = ind2sub(size(Data.StallingMatrix'),idx_to_move);
%             set(handles.edit_segno,'string',num2str(seg_no))
%             set(handles.edit_volnumber,'string',num2str(frame_no))
%             draw(hObject, eventdata, handles);
%         end
%     end
else
    set(handles.menu_validateStalls,'Checked','off')
    set(handles.uipanel_validationPanel,'Visible','off')
    delete(handles.axes4.Children)
    delete(handles.axes5.Children)
    delete(handles.axes6.Children)
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
    updateSegQualityLength(handles.axes4,str2num(handles.pixelCutofflengthValue.String)+1,'r')
    handles.NumofFilteredValue.String = num2str(sum(pixelLength <= str2num(handles.pixelCutofflengthValue.String)));
    grid on
    Data.segAnalysis.cutoff = str2num(handles.pixelCutofflengthValue.String);
    Data.segAnalysis.capNum = ind2sub(size(pixelLength),find(pixelLength <=5));
    
    if ~isfield(Data.segAnalysis,'DatasetAveInt') && ~isfield(Data.segAnalysis,'DatasetAveSNR')
        for j = 1:length(Data.seg)
            for i = 1:size(Data.I,3)
                oneRow = Data.seg(j).LRimage(:,i);
                Data.segAnalysis.AveIntensity(j,i) = mean(oneRow);
                Data.segAnalysis.AveCOV(j,i) = mean(oneRow)/std(oneRow);
            end
            Data.segAnalysis.DatasetAveInt(j,1) = mean(Data.segAnalysis.AveIntensity(j,:)/max(Data.segAnalysis.AveIntensity(j,:)));
            Data.segAnalysis.DatasetAveCOV(j,1) = mean(Data.segAnalysis.AveCOV(j,:)/max(Data.segAnalysis.AveCOV(j,:)));
        end
    end
    
function makeSegQualityAnalysis(handles)
    global Data
    delete(handles.axes5.Children)
    handles.axes5.NextPlot = 'add';
    seg_no = str2double(get(handles.edit_segno,'string'));
    LegendString = {};
    
    if handles.checkbox_IntensityDisp.Value == 1
        aveIntensity = Data.segAnalysis.AveIntensity(seg_no,:);
        IntHandle = plot(handles.axes5,1:size(Data.I,3),aveIntensity/max(aveIntensity),'-o','LineWidth',1,'MarkerSize',2);
        IntHandle.Color = 'r';
        LegendString{1,end+1} = "Normalized Ave Intensity: " + num2str(mean(aveIntensity/max(aveIntensity)));
    end
    if handles.checkbox_COVDisp.Value == 1
        COV = Data.segAnalysis.AveCOV(seg_no,:);
        COVHandle = plot(handles.axes5,1:size(Data.I,3),COV/max(COV),'-o','LineWidth',1,'MarkerSize',2);
        COVHandle.Color = 'g';
        LegendString{1,end+1} = "Normalized Ave 1/COV: "+ num2str(mean(COV/max(COV)));
    end
    handles.axes5.XAxis.Label.String = 'Slice Num';
    handles.axes5.XAxis.FontWeight = 'Bold';
    handles.axes5.XAxis.TickValues = linspace(0,350,15);
    handles.axes5.YAxis.FontSize = 10;
    handles.axes5.YAxis.FontWeight = 'Bold';
    handles.axes5.YAxis.TickValues = linspace(0,1,6);

    if ~isempty(LegendString)
        lgd = legend(handles.axes5,LegendString);
            lgd.Location = 'northoutside';
            lgd.Orientation = 'horizontal';
    end

function makeDatasetQualityAnalysis(handles)
    global Data
    delete(handles.axes6.Children)
    handles.axes6.NextPlot = 'add';
    seg_no = str2double(get(handles.edit_segno,'string'));
    
    if handles.checkbox_IntensityDisp.Value == 1
        DatasetAveInt = Data.segAnalysis.DatasetAveInt;
        DatasetIntHist = histogram(handles.axes6,DatasetAveInt,linspace(0,1,21));
            DatasetIntHist.FaceAlpha = 0.5;
            DatasetIntHist.FaceColor = 'r';
            
        updateSegQualityLength(handles.axes6,Data.segAnalysis.DatasetAveInt(seg_no,1),'r')
    end

    if handles.checkbox_COVDisp.Value == 1
        DatasetAveCOV = Data.segAnalysis.DatasetAveCOV;
        DatasetSNRHist = histogram(handles.axes6,DatasetAveCOV,linspace(0,1,21));
            DatasetSNRHist.FaceAlpha = 0.5;
            DatasetSNRHist.FaceColor = 'g';
        updateSegQualityLength(handles.axes6,Data.segAnalysis.DatasetAveCOV(seg_no,1),'g')
    end
   
        handles.axes6.XAxis.Limits = [0,1];
        handles.axes6.XAxis.TickValues = linspace(0,1,11);
        handles.axes6.XAxis.Label.String = 'Normalized Average';
        handles.axes6.XAxis.FontWeight = 'Bold';
        handles.axes6.YAxis.Label.String = 'Num of Capillary';
        handles.axes6.YAxis.FontWeight = 'Bold';

function updateSegQualityLength(axes,AveNum,Color)
    Line = xline(axes,AveNum,Color,'LineWidth',2);

function updateStallandFrame(handles, mov_dir)

global Data
seg_no = str2double(get(handles.edit_segno,'string'));
frame_no = str2double(get(handles.edit_volnumber,'string'));
SkipAllDisplay(handles);
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
        [frame_no, seg_no] = ind2sub(size(Data.StallingMatrix'),idx_to_move);
        set(handles.edit_segno,'string',num2str(seg_no))
        set(handles.edit_volnumber,'string',num2str(frame_no))
        set(handles.slider_movedata,'Value',frame_no);
        (handles);
        draw([], [], handles);
    end
end

function SkipAllDisplay(handles)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
if isfield(Data,'segAnalysis')
    if ismember(seg_no,Data.segAnalysis.capNum)
        handles.pushbutton_SkipShortSeg.Visible = 'on';
        handles.pushbutton_unSkipShortSeg.Visible = 'on';
    else
        handles.pushbutton_SkipShortSeg.Visible = 'off';
        handles.pushbutton_unSkipShortSeg.Visible = 'off';
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


% --- Executes on button press in displayCurrentSegment.
function checkbox_displayCurrentSegment_Callback(hObject, eventdata, handles)
% hObject    handle to displayCurrentSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayCurrentSegment
draw(hObject, eventdata, handles)



function pixelCutofflengthValue_Callback(hObject, eventdata, handles)
% hObject    handle to pixelCutofflengthValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelCutofflengthValue as text
%        str2double(get(hObject,'String')) returns contents of pixelCutofflengthValue as a double
makeSegLengthHistogram(handles)


% --- Executes during object creation, after setting all properties.
function pixelCutofflengthValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelCutofflengthValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_SkipShortSeg.
function pushbutton_SkipShortSeg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SkipShortSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
frame_no = str2double(get(handles.edit_volnumber,'string'));
if isfield(Data,'GTStallingMatrix')
    Data.GTStallingMatrix(seg_no,:) = 3;
    Data.ValidationFlag(seg_no,:) = 1;
end
updateStallandFrame(handles, 'Next')


% --- Executes on button press in pushbutton_unSkipShortSeg.
function pushbutton_unSkipShortSeg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_unSkipShortSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
seg_no = str2double(get(handles.edit_segno,'string'));
if isfield(Data,'GTStallingMatrix')
    Data.GTStallingMatrix(seg_no,:) = 0;
    Data.ValidationFlag(seg_no,:) = 0;
end


% --- Executes during object creation, after setting all properties.
function NumofFilteredValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumofFilteredValue (see GCBO)
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
    KeyPressed = eventdata.Key;
    if isfield(Data,'seg')
        if strcmpi(KeyPressed,"uparrow") % Previous Cap
            pushbutton_prevseg_Callback(hObject, eventdata, handles);
            disp('Previous Capillary')
        end

        if strcmpi(KeyPressed,"downarrow") % Next Cap
            pushbutton_nextseg_Callback(hObject, eventdata, handles);
            disp('Next Capillary')
        end
    end
    
    if strcmpi(KeyPressed,"leftarrow") % last image
        pushbutton_moveleft_Callback(hObject, eventdata, handles);
        disp('Previous Image')
    end

    if strcmpi(KeyPressed,"rightarrow") % next image
        pushbutton_moveright_Callback(hObject, eventdata, handles);
        disp('Next Image')
    end
    
    if strcmpi(KeyPressed,"w") % increase upper throshold of contrast
        handles.edit_MaxI.String = num2str(str2double(handles.edit_MaxI.String) + 100);
        edit_MaxI_Callback(hObject, eventdata, handles);
        disp('Increase Max Threshold')
    end
    
    if strcmpi(KeyPressed,"s") % decrease upper throshold of contrast
        handles.edit_MaxI.String = num2str(str2double(handles.edit_MaxI.String) - 100);
        edit_MaxI_Callback(hObject, eventdata, handles);
        disp('Decrease Max Threshold')
    end

    % Validation Use
    if ~strcmpi(get(handles.menu_validateStalls,'Checked'), 'off')
        if strcmpi(KeyPressed,"a") % Prev
            pushbutton_prevStallforVerification_Callback(hObject, eventdata, handles);
            disp('Previous Case for validation')
        end

        if strcmpi(KeyPressed,"d") % Next
            pushbutton_nextStallforVerification_Callback(hObject, eventdata, handles);
            disp('Next Case for validation')
        end
        
        if strcmpi(KeyPressed,"j") % Stall
            pushbutton_isaStall_Callback(hObject, eventdata, handles)
            disp('Select Stall')
        end
        
        if strcmpi(KeyPressed,"k") % Non-Stall
            pushbutton_notaStall_Callback(hObject, eventdata, handles)
            disp('Select Non-Stall')
        end
        
        if strcmpi(KeyPressed,"l") % Questionable
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
if isfield(Data,'segAnalysis')
    makeSegQualityAnalysis(handles)
    makeDatasetQualityAnalysis(handles)
else
    handles.checkbox_COVDisp.Value = 0;
    handles.checkbox_IntensityDisp.Value = 0; 
end

% --- Executes on button press in checkbox_COVDisp.
function checkbox_COVDisp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_COVDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_COVDisp
global Data
if isfield(Data,'segAnalysis')
    makeSegQualityAnalysis(handles)
    makeDatasetQualityAnalysis(handles)
else
    handles.checkbox_COVDisp.Value = 0;
    handles.checkbox_IntensityDisp.Value = 0; 
end
