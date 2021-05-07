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

% Last Modified by GUIDE v2.5 01-Mar-2021 15:33:44

% Begin initialization code - DO NOT EDIT

global TN_brain;
global TN_scalp;
global TN_roi;
global TN_brain_roi_interface;
global TN_brain_color;
global TN_scalp_color;
global TN_roi_color;
global TN_brain_around_roi_center;
global TN_roi_size;
global TN_roi_center;
global TN_roi_normal;
global TN_roi_normal_len;
global TN_scalp_normal;
global TN_roi_normal_saved;

TN_roi_normal_len=10;
TN_brain_around_roi_center=10;
TN_brain_color = [0.7 0.7 0.7];
TN_scalp_color = [0.8945 0.7578 0.5938];
TN_roi_color = [0.0 0.5 0.0];

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
       if ( strcmp(lower(varargin{i}),'roi') && i<nr_inputs) 
          if (isstruct(varargin{i+1}))
           TN_roi=varargin{i+1};
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
roi_center=mean(TN_roi.node',1);
r = get_mesh_elm_centers(TN_brain);
r(:,1)=r(:,1)-roi_center(1);
r(:,2)=r(:,2)-roi_center(2);
r(:,3)=r(:,3)-roi_center(3);
r=sqrt(sum(r.^2,2));
ind=find(r<=TN_brain_around_roi_center);
TN_brain_roi_interface=TN_brain;
TN_brain_roi_interface.face=TN_brain_roi_interface.face(:,ind);
TN_brain_roi_interface.field=zeros(1,length(ind));
TN_brain_roi_interface = clean_tri_mesh(TN_brain_roi_interface.node',TN_brain_roi_interface.face',TN_brain_roi_interface.field');

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


% --- Executes just before TargetingNavigator is made visible.
function TargetingNavigator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TargetingNavigator (see VARARGIN)

% Choose default command line output for TargetingNavigator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global TN_surf_first;
global TN_surf_second;
global TN_alpha_first;
global TN_alpha_second;
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

popup_menus3={'ROI Size','Tangential E-Field Component','Hair thickness','Coil Pitch','Coil Roll','Roil Yaw'};
handles.popupmenu3.String = popup_menus3;


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
             tmp=quiver3(TN_roi_center(1),TN_roi_center(2),TN_roi_center(3),TN_roi_normal(1)*TN_roi_normal_len,TN_roi_normal(2)*TN_roi_normal_len,TN_roi_normal(3)*TN_roi_normal_len);
             tmp.Color = 'cyan';
             tmp.LineWidth=5;
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
             tmp=quiver3(TN_roi_center(1),TN_roi_center(2),TN_roi_center(3),TN_roi_normal(1)*TN_roi_normal_len,TN_roi_normal(2)*TN_roi_normal_len,TN_roi_normal(3)*TN_roi_normal_len);
             tmp.Color = 'cyan';
             tmp.LineWidth=3;
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
