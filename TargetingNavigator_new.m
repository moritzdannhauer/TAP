function varargout = TargetingNavigator(varargin)
% TARGETINGNAVIGATOR MATLAB code for TargetingNavigator.fig
%      TARGETINGNAVIGATOR, by itself, creates a new TARGETINGNAVIGATOR or raises the existing
%      singleton*.
%
%      H = TARGETINGNAVIGATOR returns the handle to a new TARGETINGNAVIGATOR or the handle to
%      the existing singleton*.
%
%      TARGETINGNAVIGATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TARGETINGNAVIGATOR.M with the given input arguments.
%
%      TARGETINGNAVIGATOR('Property','Value',...) creates a new TARGETINGNAVIGATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TargetingNavigator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TargetingNavigator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TargetingNavigator

% Last Modified by GUIDE v2.5 29-Apr-2022 03:23:33

% Begin initialization code - DO NOT EDIT
global tmp;
global TN_brain;
global TN_scalp;
global TN_brain_tetcenter;
global TN_brain_roi_interface;
global TN_brain_color;
global TN_scalp_color;
global TN_roi_color;
global TN_brain_around_roi_center;
global TN_roi;
global TN_roi_size;
global TN_roi_center;
global TN_roi_normal;
global TN_roi_normal_len;
global TN_scalp_normal;
global TN_roi_normal_saved;
global TN_selected_stimulator;
global TN_selected_coil;
global TN_chosen_stimulator_coil_setup;

global r30_b65_a;
global r30_b65_b;
global r30_b70_a;
global r30_b70_b;
global r30_b80_a;
global r30_b80_b;
global r30_b91_a;
global r30_b91_b;

global x100_b65_a;
global x100_b65_b;
global x100_b70_a;
global x100_b70_b;
global x100_b80_a;
global x100_b80_b;
global x100_b91_a;
global x100_b91_b;

tmp=NaN;

r30_b65_a = 1.5743;
r30_b65_b = -2.533;
r30_b70_a = 1.4627;
r30_b70_b = -1.7228;
r30_b80_a = 1.5591;
r30_b80_b = -2.1386;
r30_b91_a = 1.9308;
r30_b91_b = -2.4018;

x100_b65_a = 1.5779;
x100_b65_b = -2.856;
x100_b70_a = 1.493;
x100_b70_b = -2.0316;
x100_b80_a = 1.5442;
x100_b80_b = -2.1018;
x100_b91_a = 1.9204;
x100_b91_b = -2.9351;

TN_roi_normal_len=10;
TN_brain_around_roi_center=10;
TN_brain_color = [0.7 0.7 0.7];
TN_scalp_color = [0.8945 0.7578 0.5938];
TN_roi_color = [0.0 0.5 0.0];

TN_chosen_stimulator_coil_setup=1;
i=1; 
nr_inputs=length(varargin);
while nr_inputs>0
    if (i>=nr_inputs+1)
      break;  
    end
    if isstr(varargin{i})
       if ( strcmp(lower(varargin{i}),'brain') && i<nr_inputs) 
          if (isstruct(varargin{i+1}))
           TN_brain=varargin{i+1};
           i=i+1;
          end
       end 
       if ( strcmp(lower(varargin{i}),'scalp') && i<nr_inputs) 
          if (isstruct(varargin{i+1}))
           TN_scalp=varargin{i+1};
           i=i+1;
          end
       end 
       if ( strcmp(lower(varargin{i}),'roi_center') && i<nr_inputs) 
           TN_roi_center=varargin{i+1};
           i=i+1;
       end 
       if ( strcmp(lower(varargin{i}),'brain_tetcenter') && i<nr_inputs) 
          if (isstruct(varargin{i+1}))
           TN_brain_tetcenter=varargin{i+1};
           i=i+1;
          end
       end    
       if ( strcmp(lower(varargin{i}),'roi_size') && i<nr_inputs) 
          if (isvector(varargin{i+1}))
           TN_roi_size=varargin{i+1};
           i=i+1;
          end
       end
       if ( strcmp(lower(varargin{i}),'roi_center') && i<nr_inputs)
          if (isvector(varargin{i+1}))
           TN_roi_center=varargin{i+1};
           i=i+1;
          end
       end 
       
       if ( strcmp(lower(varargin{i}),'roi_normal') && i<nr_inputs)
          if (isvector(varargin{i+1}))
           TN_roi_normal=-1*varargin{i+1}; %Scale it to be visable, outward-pointing
           TN_roi_normal=TN_roi_normal/norm(TN_roi_normal);
           TN_roi_normal_saved=TN_roi_normal;
           i=i+1;
          end
       end 
       if ( strcmp(lower(varargin{i}),'scalp_normal') && i<nr_inputs)
          if (isvector(varargin{i+1}))
           TN_scalp_normal=varargin{i+1}; %Scale it to be visable, outward-pointing
           TN_scalp_normal=TN_scalp_normal/norm(TN_scalp_normal);
           i=i+1;
          end
       end 
       
    end
    i=i+1;
end

%brain cut with ROI
r = TN_brain_tetcenter.field';
r(:,1)=r(:,1)-TN_roi_center(1);
r(:,2)=r(:,2)-TN_roi_center(2);
r(:,3)=r(:,3)-TN_roi_center(3);
r=sqrt(sum(r.^2,2));
ind=find(r<=TN_brain_around_roi_center);
TN_brain_roi_interface=TN_brain;
TN_brain_roi_interface.face=TN_brain_roi_interface.face(:,ind);
TN_brain_roi_interface.field=zeros(1,length(ind));
TN_brain_roi_interface = clean_tri_mesh(TN_brain_roi_interface.node',TN_brain_roi_interface.face',TN_brain_roi_interface.field');

ind=find(r<=TN_roi_size*2);
TN_roi.node=TN_brain_tetcenter.node';
TN_roi.face=TN_brain_tetcenter.face(:,ind);
TN_roi.field=zeros(1,length(ind));
TN_roi = clean_tri_mesh(TN_roi.node',TN_roi.face',TN_roi.field');

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TargetingNavigator_OpeningFcn, ...
                   'gui_OutputFcn',  @TargetingNavigator_OutputFcn, ...
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

function TargetingNavigator_PlotSurf(varargin)

if (mod(length(varargin),3)==0)
  cla reset;  
 for i=1:3:length(varargin)   
  tmp=varargin{i};
  points=tmp.node';
  faces=tmp.face';
  face_color = varargin{i+1};
  face_alpha = varargin{i+2};
  cmap = gray(256);
  edge_color = 'none';
  culling = false;
  backfacelighting = 'lit';
  handles = patch('vertices',points,'faces',faces,'facecolor',face_color,'edgecolor', ... 
            edge_color,'backfacelighting',backfacelighting,'FaceAlpha',face_alpha);
  lighting phong;
  material([ 0.00 0.50 0.20 2.00 1.00 ]);
 end
  hl(1) = camlight(0,40,'infinite');
  hl(2) = camlight(180,40,'infinite');
  hl(3) = camlight(0,-90,'infinite');
  hl(4) = camlight(90,0,'infinite');
  hl(4) = camlight(-90,0,'infinite');

 camproj('perspective');
 axis square;
 axis off;
 axis tight;
 axis equal;
end

function populate_tms_intensity(selection)
 global handles;
 global TN_chosen_stimulator_coil_setup;
 global r30_b65_a;
 global r30_b65_b;
 global r30_b70_a;
 global r30_b70_b;
 global r30_b80_a;
 global r30_b80_b;
 global r30_b91_a;
 global r30_b91_b;

 global x100_b65_a;
 global x100_b65_b;
 global x100_b70_a;
 global x100_b70_b;
 global x100_b80_a;
 global x100_b80_b;
 global x100_b91_a;
 global x100_b91_b;

 global MSO;
 global didt;
 
 TN_chosen_stimulator_coil_setup = selection;
 
 switch selection
      case 1 
           a=r30_b65_a;
           b=r30_b65_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b; 
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));         
      case 2
           a=r30_b70_a;
           b=r30_b70_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b; 
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));  
      case 3
           a=r30_b80_a;
           b=r30_b80_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b;   
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 4
           a=r30_b91_a;
           b=r30_b91_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b; 
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 5
           a=x100_b65_a;
           b=x100_b65_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b;  
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 6
           a=x100_b70_a;
           b=x100_b70_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b;  
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 7
           a=x100_b80_a;
           b=x100_b80_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b;  
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 8
           a=x100_b91_a;
           b=x100_b91_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           didt=MSO*a+b;  
           if (didt<0)
             didt=0;  
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
 end

function populate_tms_intensity_didt_changed(selection)
 global handles;
 global TN_chosen_stimulator_coil_setup;
 global r30_b65_a;
 global r30_b65_b;
 global r30_b70_a;
 global r30_b70_b;
 global r30_b80_a;
 global r30_b80_b;
 global r30_b91_a;
 global r30_b91_b;

 global x100_b65_a;
 global x100_b65_b;
 global x100_b70_a;
 global x100_b70_b;
 global x100_b80_a;
 global x100_b80_b;
 global x100_b91_a;
 global x100_b91_b;

 global MSO;
 global didt;
 global tmp;
 TN_chosen_stimulator_coil_setup = selection;
 
 switch selection
      case 1 
           a=r30_b65_a;
           b=r30_b65_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a); 
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               didt=tmp;
               MSO=val;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));
 
      case 2
           a=r30_b70_a;
           b=r30_b70_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a); 
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 3
           a=r30_b80_a;
           b=r30_b80_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a);    
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));
      case 4
           a=r30_b91_a;
           b=r30_b91_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a);  
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));
      case 5
           a=x100_b65_a;
           b=x100_b65_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a); 
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO));
      case 6
           a=x100_b70_a;
           b=x100_b70_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a); 
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 7
           a=x100_b80_a;
           b=x100_b80_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a);  
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
      case 8
           a=x100_b91_a;
           b=x100_b91_b;
           set(handles.edit17,'String',num2str(a));
           set(handles.edit18,'String',num2str(b));
           val=round((didt-b)/a); 
           if (~isnan(tmp))
             val=round((tmp-b)/a);
           end
           if (val<=100 && val>=0)
               MSO=val;
               didt=tmp;
           end
           set(handles.edit15,'String',num2str(didt));
           set(handles.edit16,'String',num2str(MSO)); 
 end 
 tmp=NaN;
 
function TargetingNavigator_DeleteFcn(hObject, eventdata, handles, varargin)
  
disp('done');

% --- Executes just before TargetingNavigator is made visible.
function TargetingNavigator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TargetingNavigator (see VARARGIN)
global handles;
% Choose default command line output for TargetingNavigator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global TN_surf_first;
global TN_surf_second;
global TN_alpha_first;
global TN_alpha_second;
global TN_selected_stimulator;
global TN_selected_coil;

global x100_b65_a;
global x100_b65_b;

global popup_menus9;
global popup_menus10;


TN_surf_first='Brain Surface';
TN_surf_second='None';
TN_alpha_first=0.5;
TN_alpha_second=1.0;

%Setup menus
popup_menus={'Brain Surface','Brain/ROI Interface','Scalp Surface','ROI','None'};
handles.popupmenu2.String = popup_menus;
handles.popupmenu4.String = popup_menus;
handles.edit6.String = 0.5;
handles.edit7.String = 1.0;
handles.popupmenu4.Value=length(popup_menus);

%setup more menus
popup_menus2={'ROI Location','Coil Center','RPA','LPA','Nasion'};
handles.popupmenu1.String = popup_menus2;

popup_menus3={'ROI Size (Radius in mm)','Tangential E-Field Component','Hair thickness (in mm)','Coil Pitch','Coil Roll','Roil Yaw'};
handles.popupmenu3.String = popup_menus3;

popup_menus4={'Median','75% max','90% max','99% max','E20','E50','E100'};
handles.popupmenu8.String = popup_menus4;

popup_menus11={'Magnitude','E-field into wall '};
handles.popupmenu11.String = popup_menus11;

popup_menus9={'MagVenture MagPro R30','MagVenture MagPro X100'};
handles.popupmenu9.String = popup_menus9;
TN_selected_stimulator = popup_menus9(1);

popup_menus10={'A/P Cool B65','MC B70','D-B80','MRI-B91 AIR COOLED'};
handles.popupmenu10.String = popup_menus10;
TN_selected_coil = popup_menus10(1);

global MSO;
global didt;
MSO=50;
didt=MSO*x100_b65_a+x100_b65_b;
handles.edit16.String = num2str(MSO);
handles.edit15.String = num2str(didt);

handles.edit17.String = num2str(x100_b65_a);
handles.edit18.String = num2str(x100_b65_b);

function [Output1, Output2] = TargetingNavigator_getSurfAndColor(TN_surf)

global TN_brain_color;
global TN_brain;
global TN_scalp_color;
global TN_scalp;
global TN_brain_roi_interface_color;
global TN_brain_roi_interface;

Output1=[];
Output2=[];

if (~isempty(TN_surf))
  switch TN_surf
      case 'None' 
         Output1=[];
         Output2=[];
      case 'Brain Surface'
         if (exist('TN_brain','var') && exist('TN_brain_color','var'))
             if (isstruct(TN_brain) && length(TN_brain_color)==3)
                Output1=TN_brain;
                Output2=TN_brain_color;
             end
         end
      case 'Scalp Surface'
         if (exist('TN_scalp','var') && exist('TN_scalp_color','var'))
              if (isstruct(TN_scalp) && length(TN_scalp_color)==3)
                Output1=TN_scalp;
                Output2=TN_scalp_color;  
             end
         end
          
      case 'ROI'
         if (exist('TN_roi','var') && exist('TN_roi_color','var'))
              if (isstruct(TN_roi) && length(TN_roi_color)==3)
                Output1=TN_roi;
                Output2=TN_roi_color;  
             end
         end          
      case 'Brain/ROI Interface'
         if (exist('TN_brain_roi_interface','var') && exist('TN_brain_color','var'))
              if (isstruct(TN_roi) && length(TN_roi_color)==3)
                Output1=TN_brain_roi_interface;
                Output2=TN_brain_color;  
             end
         end
  end
end

% --- Outputs from this function are returned to the command line.
function varargout = TargetingNavigator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global TN_surf_first;
global TN_surf_second;
global TN_alpha_first;
global TN_alpha_second;

plotted=false;
if (exist('TN_surf_first','var'))
  [first, first_color] = TargetingNavigator_getSurfAndColor(TN_surf_first);
end

if (exist('TN_surf_second','var'))
  [second, second_color] = TargetingNavigator_getSurfAndColor(TN_surf_second); 
end

if (~isempty(first) && ~isempty(first_color) && ~isempty(second) && ~isempty(second_color))
   TargetingNavigator_PlotSurf(first,first_color,TN_alpha_first,second,second_color,TN_alpha_second);
   plotted=true;
else
   if (~isempty(first) && ~isempty(first_color))
      TargetingNavigator_PlotSurf(first,first_color,TN_alpha_first);
      plotted=true;
   end
   if (~isempty(second) && ~isempty(second_color))
      TargetingNavigator_PlotSurf(second,second_color,TN_alpha_second);
      plotted=true;
   end
end

if (~plotted)
    disp('not plotted');
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(~, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TargetingNavigator_PlotFirst()

global TN_brain;
global TN_brain_color;
global TN_scalp;
global TN_scalp_color;
global TN_surf_first;
global TN_surf_second;
global TN_roi;
global TN_roi_color;
global TN_alpha_first;
global TN_alpha_second;
global TN_brain_roi_interface;
global TN_roi_size;
global TN_roi_center;
global TN_roi_normal;
global TN_roi_normal_len;

switch TN_surf_first
    case 'None'
               switch TN_surf_second
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end  
                case 'Brain Surface'
                    if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_second); 
                    end
                case 'Scalp Surface'
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_second); 
                    end
                   case 'None'
                       cla reset; 
                otherwise
                    
             end
    case 'Brain Surface'
        if (isstruct(TN_brain) && ~isempty(TN_brain_color))
         if(strcmp(TN_surf_second,'None'))   
          TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first);
         else
            switch TN_surf_second
                case 'Brain Surface'

                case 'Scalp Surface'
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first,TN_roi,TN_roi_color,TN_alpha_second); 
                    end
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_second,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end  
                otherwise
                    
             end
         end
        end
    case 'Scalp Surface'
        if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
         if(strcmp(TN_surf_second,'None'))   
          TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first);
         else
            switch TN_surf_second
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_second,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end  
                case 'Brain Surface'
                    if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_second,TN_brain,TN_brain_color,TN_alpha_second); 
                    end
                case 'Scalp Surface'

                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first,TN_roi,TN_roi_color,TN_alpha_second); 
                    end
                otherwise
                    
             end
         end
        end
    case 'ROI'
        if(strcmp(TN_surf_second,'None'))   
          TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first);
        else
          switch TN_surf_second
             case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end  
              case 'ROI'
                  
              case 'Scalp Surface'    
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
              case 'Brain Surface'    
                     if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain,TN_brain_color,TN_alpha_second); 
                    end
              otherwise
                  
          end
        end
        
        case 'Brain/ROI Interface'
        if(strcmp(TN_surf_second,'None'))   
          TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first);
        else
          switch TN_surf_second
             case 'Brain/ROI Interface'
 
              case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end
              case 'Scalp Surface'    
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
              case 'Brain Surface'    
                     if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain,TN_brain_color,TN_alpha_second); 
                    end
              otherwise
                  
          end
        end  
    otherwise
        disp('??')  
end 
if(strcmp(TN_surf_second,'ROI') || strcmp(TN_surf_first,'ROI'))   
     if (exist('TN_roi_center','var') && exist('TN_roi_normal','var'))
         tmp=size(TN_roi_center);
         if (min(tmp)==1 && max(tmp)==3)
          tmp=size(TN_roi_normal);   
           if (min(tmp)==1 && max(tmp)==3)
             hold on;
             tmp=quiver3(TN_roi_center(1)+TN_roi_normal_len*TN_roi_normal(1),TN_roi_center(2)+TN_roi_normal_len*TN_roi_normal(2),TN_roi_center(3)+TN_roi_normal_len*TN_roi_normal(3),-TN_roi_normal(1)*TN_roi_normal_len,-TN_roi_normal(2)*TN_roi_normal_len,-TN_roi_normal(3)*TN_roi_normal_len);
             tmp.Color = 'cyan';
             tmp.LineWidth=5;
             tmp.MaxHeadSize=3;
           end
         end
     end
end  

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)

% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

global TN_surf_first;

contents = cellstr(get(hObject,'String'));
TN_surf_first=contents{get(hObject,'Value')};

TargetingNavigator_PlotFirst();


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global TN_rotate_that;
global TN_scalp_normal;
global TN_roi_normal;
global TN_roi_normal_saved;

keep_old_value=false;
tmp=str2double(handles.edit5.String);

if (isnumeric(tmp))
    if (tmp>=-180 && tmp<=180)
        TN_rotate_that=tmp;
    else
        keep_old_value=true;
    end
end

if (~keep_old_value && strcmp(char(handles.popupmenu3.String(handles.popupmenu3.Value)),'Tangential E-Field Component'))
 TN_roi_normal = rotVecAroundArbAxis(TN_roi_normal_saved,TN_scalp_normal,tmp);
 cla reset;
 TargetingNavigator_PlotFirst();
 TargetingNavigator_PlotSecond();
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1



function edit6_Callback(hObject, eventdata, handles)

global TN_alpha_first;

keep_old_value=false;
tmp=str2double(handles.edit6.String);
if (isnumeric(tmp))
    if (tmp>=0.02 && tmp<=1.0)
        TN_alpha_first=tmp;
    else
        keep_old_value=true;
    end
else
    keep_old_value=true;
end

if (keep_old_value)
   handles.edit6.String=num2str(TN_alpha_first);
end

if (~keep_old_value)
  TargetingNavigator_PlotFirst();
end



% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TargetingNavigator_PlotSecond()
global TN_brain;
global TN_brain_color;
global TN_scalp;
global TN_scalp_color;
global TN_surf_first;
global TN_surf_second;
global TN_roi;
global TN_roi_color;
global TN_alpha_first;
global TN_alpha_second;
global TN_brain_roi_interface;
global TN_roi_size;
global TN_roi_center;
global TN_roi_normal;
global TN_roi_normal_len;

switch TN_surf_second
    case 'Brain/ROI Interface'
              switch TN_surf_first
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first); 
                    end 
                case 'Brain Surface'
                    if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end
                case 'Scalp Surface'
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end
                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain_roi_interface,TN_brain_color,TN_alpha_second); 
                    end
                case 'None'
                    cla reset; 
                otherwise
                    
              end
        
    case 'None'
        switch TN_surf_first
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first); 
                    end 
                case 'Brain Surface'
                    if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first); 
                    end
                case 'Scalp Surface'
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first); 
                    end
                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first); 
                    end
                case 'None'
                    cla reset; 
                otherwise
                    
             end
    case 'Brain Surface'
        if (isstruct(TN_brain) && ~isempty(TN_brain_color))
         if(strcmp(TN_surf_first,'None'))   
          TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_second);
         else
            switch TN_surf_first
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first,TN_brain,TN_brain_color,TN_alpha_second); 
                    end 
                case 'Brain Surface'

                case 'Scalp Surface'
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first,TN_brain,TN_brain_color,TN_alpha_second); 
                    end
                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_brain,TN_brain_color,TN_alpha_second); 
                    end
                otherwise
                    
             end
         end
        end
    case 'Scalp Surface'
        if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
         if(strcmp(TN_surf_first,'None'))   
          TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_second);
         else
            switch TN_surf_first
                case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end 
                case 'Brain Surface'
                    if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
                case 'Scalp Surface'

                case 'ROI'
                    if (isstruct(TN_roi) && ~isempty(TN_roi_color))
                       TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_first,TN_scalp,TN_scalp_color,TN_alpha_second); 
                    end
                otherwise
                    
             end
         end
        end
    case 'ROI'
        if(strcmp(TN_surf_first,'None'))   
          TargetingNavigator_PlotSurf(TN_roi,TN_roi_color,TN_alpha_second);
        else
          switch TN_surf_first
              case 'Brain/ROI Interface'
                    if (isstruct(TN_brain_roi_interface) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain_roi_interface,TN_brain_color,TN_alpha_first,TN_roi,TN_roi_color,TN_alpha_second); 
                    end 
              case 'ROI'
                  
              case 'Scalp Surface'    
                    if (isstruct(TN_scalp) && ~isempty(TN_scalp_color))
                       TargetingNavigator_PlotSurf(TN_scalp,TN_scalp_color,TN_alpha_first,TN_roi,TN_roi_color,TN_alpha_second); 
                    end
              case 'Brain Surface'    
                     if (isstruct(TN_brain) && ~isempty(TN_brain_color))
                       TargetingNavigator_PlotSurf(TN_brain,TN_brain_color,TN_alpha_first,TN_roi,TN_roi_color,TN_alpha_second); 
                    end
              otherwise
                  
          end
        end
    otherwise
        disp('??')
end 
if(strcmp(TN_surf_second,'ROI') || strcmp(TN_surf_first,'ROI'))   
     if (exist('TN_roi_center','var') && exist('TN_roi_normal','var'))
         tmp=size(TN_roi_center);
         if (min(tmp)==1 && max(tmp)==3)
          tmp=size(TN_roi_normal);   
           if (min(tmp)==1 && max(tmp)==3)
             hold on;
             tmp=quiver3(TN_roi_center(1)+TN_roi_normal_len*TN_roi_normal(1),TN_roi_center(2)+TN_roi_normal_len*TN_roi_normal(2),TN_roi_center(3)+TN_roi_normal_len*TN_roi_normal(3),-TN_roi_normal(1)*TN_roi_normal_len,-TN_roi_normal(2)*TN_roi_normal_len,-TN_roi_normal(3)*TN_roi_normal_len);
             tmp.Color = 'cyan';
             tmp.LineWidth=5;
             tmp.MaxHeadSize=3;
           end
         end
     end
end

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

global TN_surf_second;

contents = cellstr(get(hObject,'String'));
TN_surf_second=contents{get(hObject,'Value')};

TargetingNavigator_PlotSecond();


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

global TN_alpha_second;

keep_old_value=false;
tmp=str2double(handles.edit7.String);
if (isnumeric(tmp))
    if (tmp>=0.02 && tmp<=1.0)
        TN_alpha_second=tmp;
    else
        keep_old_value=true;
    end
else
    keep_old_value=true;
end

if (keep_old_value)
   handles.edit7.String=num2str(TN_alpha_second);
end

if (~keep_old_value)
  TargetingNavigator_PlotSecond();
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d on;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 zoom on;



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

global TN_cut_brain_radius;
TN_cut_brain_radius=str2double(get(hObject,'String'))

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
global didt;
global MSO;
global TN_chosen_stimulator_coil_setup;
global tmp;
val = str2double(get(hObject,'String'));
if ( val == 0)
  MSO=0; 
  didt=0;
  set(handles.edit15,'String',num2str(didt));
  set(handles.edit16,'String',num2str(MSO)); 
  return ;
end
if (val<0)
   set(handles.edit15,'String',num2str(didt));
   set(handles.edit16,'String',num2str(MSO)); 
else
   tmp=val;
   populate_tms_intensity_didt_changed(TN_chosen_stimulator_coil_setup);
end 
 

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
global MSO;
global didt;
global TN_chosen_stimulator_coil_setup;
val = round(str2double(get(hObject,'String')));

if (val>0 && val<=100)
 MSO = round(str2double(get(hObject,'String')));
 populate_tms_intensity(TN_chosen_stimulator_coil_setup);
 return ;
end
    
if ( val == 0 )
 didt=0;  
 MSO=0;
end

 set(handles.edit15,'String',num2str(didt));
 set(handles.edit16,'String',num2str(MSO)); 




% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
global TN_Efield_metric;

TN_Efield_metric=get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
global TN_Desired_ROI_Efield_Value;
TN_Desired_ROI_Efield_Value = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9
global TN_selected_stimulator;
global TN_selected_coil;
global popup_menus9;
global popup_menus10;

factor = get(hObject,'Value');
TN_selected_stimulator=popup_menus9(factor);
count = 1;

for i=1:length(popup_menus10)
    if(strcmp(popup_menus10(i),TN_selected_coil))
      count=i;  
      break;
    end
end

factor=(length(popup_menus10)*(factor-1))+count;
populate_tms_intensity(factor);

% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10
global TN_selected_stimulator;
global TN_selected_coil;
global popup_menus9;
global popup_menus10;

factor = get(hObject,'Value');
TN_selected_coil = popup_menus10(factor);
count = 1;

for i=1:length(popup_menus9)
    if(strcmp(popup_menus9(i),TN_selected_stimulator))
      count=i;  
      break;
    end
end

factor=(length(popup_menus10)*(count-1))+factor;

populate_tms_intensity(factor);

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11
global TN_EfieldScaling;

TN_EfieldScaling=get(hObject,'Value');

 
% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function uibuttongroup4_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global abc;
abc=22;
disp('aus');
