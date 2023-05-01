clear all;

global root_folder;

if (ispc)
    sep='\';
    not_sep='/';
    rep_space = ' ';
elseif (ismac || isunix)
    sep='/';
    not_sep='\';
    rep_space = '\ ';
end
space = ' ';

%%%%% CHANGE TAP PARAMETER HERE: Begin %%%%%
MSO=60; %this maximum stimulator output could be, for example, resting motor threshold
Scale_Efield_to_Value=81.8; % NaN = scales to MSO, otherwise to the given value
NO_USE_THIS_MAX_DIDT=NaN;
Efield_display_decimal_places=2;
Suppress_Any_Konsole_Input_or_Graphical_Output=1;
Stokes_factor=2.7;
Neuronaviagtion_system='Brainsight'; %'Localite' 'ANT Visor2'
DTI=0;
StartTargetingNavigator = 0;
fsl_path='/usr/local/fsl';
simnibs_folder='/Users/m/Applications/SimNIBS-3.2/';
root_folder=[pwd sep];
freesurfer_matlab_folder='/Applications/freesurfer/matlab/';
coil_models_are_here='/simnibs_env/lib/python3.7/site-packages/simnibs/ccd-files/';
subj='ernie';
subjects_folder='subjects';
target_name='M1';
mask='M1_mask_TMS_target.nii.gz';
outputfolder='M1';
scirunfolder='M1_scirun';
Original_MRI='org/ernie_T1.nii.gz';
SimNIBS_MRI=[subj '_T1fs_conform.nii.gz'];
prefix_for_voxelized_nii_files='res';
dti_file = ['d2c_' subj '/dti_results_T1space/DTI_conf_tensor.nii.gz'];
%change separators and spaces
fsl_path=strrep(fsl_path,not_sep,sep);
fsl_path=strrep(fsl_path,space,rep_space);
simnibs_folder=strrep(simnibs_folder,not_sep,sep);
simnibs_folder=strrep(simnibs_folder,space,rep_space);
root_folder=strrep(root_folder,not_sep,sep);
root_folder=strrep(root_folder,space,rep_space);
coil_models_are_here=strrep(coil_models_are_here,not_sep,sep);
coil_models_are_here=strrep(coil_models_are_here,space,rep_space);
Original_MRI=strrep(Original_MRI,not_sep,sep);
Original_MRI=strrep(Original_MRI,space,rep_space);
addpath(char([simnibs_folder sep 'matlab']));
addpath(char([root_folder sep 'matlab']));
addpath(char(freesurfer_matlab_folder));
addpath(char(fsl_path));
tms_opt = opt_struct('TMSoptimize');

Optimize_ROI_Efield_Magnitude = 0; %0 means TMS coil placement is optimized for a particular ROI-E-field direction (i.e., perpendicular to sulcal wall)
hairthicknesses=0.5;%1:8;%0.5:0.5:3.0;%:0.5:7.5; %in mm
TMS_induced_Efield_flows_into_sulcalwall=1; %1=into wall, most effcicent for neuronal recruitment
scalp_search_radius=1; %25mm radius around scalp-prpjected point
scalp_coil_center_search_grid=1; %1mm Manhattan distance
coil_rotation_discretization = 90; %1 degree
tms_opt.search_angle = 180;
coil_name='Magstim_70mm_Fig8.ccd'; %Note, SimNIBS coil models are made for single pulses meaning the induced E-field points towards TMS coil handle
tms_opt.target_size=3; %radius [mm] of ROI
Automatric_Flip_current_direction= 1; %so that TMS coil cable does not dangle in subjects face
Use_original_mask_for_Efield_eval = 0; %otherwise it uses SimNIBS Target with ...
Voxelize_Efield_mask_interpolation_thres=0.1; % ... this threshold

if (DTI==1)
  tms_opt.fname_tensor=[root_folder sep subjects_folder sep subj sep dti_file];
  tms_opt.anisotropy_type='vn';
  tms_opt.aniso_maxratio=10;
  tms_opt.aniso_maxcond=2;
  if (~exist([root_folder sep subjects_folder sep subj sep dti_file]))
     disp(['[TAP] The DTI dataset ' tms_opt.fname_tensor ' could not be found. Check Path. Aborted.' ]);
     return;
  end
else
  tms_opt.anisotropy_type='scalar';
end
%%%%% CHANGE TAP PARAMETER HERE: End %%%%%

if (~strcmp(Neuronaviagtion_system,'Brainsight'))
    disp('[TAP] Currently, TAP implements support only for the Brainsight Neuronavigation system.');
    return ;
end

if (~(strcmp(coil_name(end-6:end),'.nii.gz') || strcmp(coil_name(end-3:end),'.nii') || strcmp(coil_name(end-3:end),'.ccd')))
     disp('[TAP] You need to specify either a coil file extension as *.nii.gz or *.ccd in coil_name. Aborted.');  
     return ;
end

max_didt=1;
tms_opt.method = 'ADM';
if (~strcmp(coil_name(end-3:end),'.ccd'))
  tms_opt.method = 'direct';
  disp('[TAP] Using *.nii.gz coil files will take a VERY LONG TIME! You can use *.ccd coil files with ADM which is very fast instead.'); 
else
    if (Optimize_ROI_Efield_Magnitude==0)
       tms_opt.search_angle = 360;
    else
       tms_opt.search_angle = 180;
    end
    [canonical_coil, max_didt_from_coil]=load_coil([simnibs_folder coil_models_are_here coil_name]);%load coil definition so a transfermatrix to the optimal setup can be computed
    if (isnan(max_didt_from_coil))
        max_didt=1;
    else
     max_didt=max_didt_from_coil;
    end
    if (~isnan(NO_USE_THIS_MAX_DIDT))
     max_didt=NO_USE_THIS_MAX_DIDT;
    end
end

if (max_didt==1)
  if (strfind(coil_name,'MagVenture_MC_B70'))
    max_didt=155.3;
    disp('[TAP] No di/dt_max found in ccd file. Use approximate value from Drakaki et al. 2022 but E-field intensity scaling could be off.');
  elseif (strfind(coil_name,'Magstim_70mm_Fig8'))
    max_didt=114.7;
    disp('[TAP] No di/dt_max found in coil model. Use approximate value from Drakaki et al. 2022 but E-field intensity scaling could be off.');
  else
    disp('[TAP] No di/dt_max found in coil model. The E-field intensity scaling will be off, please specify NO_USE_THIS_MAX_DIDT or use a new coil model (Drakaki et al., 2022, SimNIBS4 repository) otherwise.');
  end
 if (strfind(tms_opt.method,'direct'))
   disp('[TAP] Cannot find information on intensity scaling (di/dt->%MSO) for chosen coil model. Check Dannhauer et al., JNE, 2022 on how to do that to set NO_USE_THIS_MAX_DIDT variable.');
 end
 disp('[TAP] new *.ccd coil models available at https://github.com/simnibs/simnibs/tree/master/simnibs/resources/coil_models/Drakaki_BrainStim_2022');
 end
 disp(['[TAP] di/dt_max = ' num2str(max_didt) ' A/' sprintf('\x3BC')  's']);

if (exist('hairthicknesses','var'))
if (~isempty(hairthicknesses) && length(hairthicknesses)==1)
 if (exist([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',hairthicknesses(1))]))
     if (Suppress_Any_Konsole_Input_or_Graphical_Output==1)
        delete_prev_folders='y';
    else
       delete_prev_folders=input('[TAP] INPUT NEEDED: Do you want to delete those folders and restart (y/n, default n)? ', 's');
    end
    if (~isempty(delete_prev_folders) && strcmpi(delete_prev_folders,'y'))
      [SUCCESS,MESSAGE,MESSAGEID]=rmdir([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',hairthicknesses(1))],'s');
      [SUCCESS,MESSAGE,MESSAGEID]=rmdir([subjects_folder sep subj sep target_name  '_scirun'],'s');
      [SUCCESS,MESSAGE,MESSAGEID]=rmdir([subjects_folder sep subj sep target_name ],'s');
    else
       disp('[TAP] TAP already ran on your single defined hairthickness. Abort.'); 
       return; 
    end
 end
end
else
    disp('[TAP] The variable hairthickness is not defined. Abort.');  
    return;
end

[brain, brain_tet_mesh, brain_surf_tetcenter, scalp, head_model_mesh, brain_vert_neigh]=get_meshSurfs(subjects_folder, subj, sep);
if ( isempty(brain) || isempty(brain_tet_mesh) || isempty(brain_surf_tetcenter) || isempty(scalp) || isempty(head_model_mesh) || isempty(brain_vert_neigh) ) 
  return ;
end
[coord, ROI_size, which_pipeline, ROI_mask_vol_file, ROI_voxel_centers] = get_target(subjects_folder,subj,mask,SimNIBS_MRI, sep, brain, Original_MRI, fsl_path);

if (ROI_size<=0 || length(coord)~=3 || ~(size(ROI_voxel_centers,1)>0 && size(ROI_voxel_centers,2)==3))
  return ;  
end
r = get_mesh_elm_centers(brain_tet_mesh);
[mindis index dis] = min_distance(r, coord);
if (mindis<ROI_size)
    disp(['[TAP] INFO only - Automatic estimation of ROI radius size determined to be: ' num2str(ROI_size) ' mm, however, the user-specified value ' num2str(tms_opt.target_size) ' is used. ']);
    disp(['[TAP] INFO only - Automatic estimation of ROI radius size: it is better to use a larger ROI because for too small values no GM tetrahedron may be found. ']);
else
    disp(['[TAP] INFO only - Automatic estimation of ROI radius is ' num2str(ROI_size) ' mm, smaller than distance to next GM voxel ' num2str(mindis) ' mm so better use ' num2str(tms_opt.target_size) ' mm you specified instead !?']);
    if (Suppress_Any_Konsole_Input_or_Graphical_Output==0)
     ROI_size=input('[TAP] Please enter ROI size (equals radius of a spherical ROI) [mm]: ', 's');
    else
     ROI_size=tms_opt.target_size;  
    end
end

min_ROI_voxel_centers=min(ROI_voxel_centers);
max_ROI_voxel_centers=max(ROI_voxel_centers);
diff=norm([min_ROI_voxel_centers-max_ROI_voxel_centers]);
tmp=r;
tmp(:,1)=tmp(:,1)-coord(1);
tmp(:,2)=tmp(:,2)-coord(2);
tmp(:,3)=tmp(:,3)-coord(3);
tmp=sqrt(sum(tmp.^2,2));
index=find(tmp<diff);

%[ind] = min_distance_thres(r(index,:), ROI_voxel_centers, 0.5);
%ROI_tet_number=brain_tet_mesh.field(index(unique(ind)));

if (isnan(coord))
    disp('[TAP] Target coordinate could be determined. Abort.');
    return;
end

[target, scalp_normal] = get_target_normal(subjects_folder, sep, subj, coord, brain, scalp, head_model_mesh, which_pipeline, brain_vert_neigh);

if (~exist([subjects_folder sep subj sep scirunfolder])) 
  mkdir([subjects_folder sep subj sep scirunfolder]);
end

if (~exist([subj '/' outputfolder ])) 
  mkdir([subjects_folder sep subj sep outputfolder ]);
end

stlwrite([subjects_folder sep subj sep outputfolder sep subj '_LoadAsNIFTI_brain.stl'],brain.face',brain.node');
stlwrite([subjects_folder sep subj sep outputfolder sep subj '_LoadAsNIFTI_scalp.stl'],scalp.face',scalp.node'); 

if (which_pipeline==3)
 tms_opt.fnamehead = [subjects_folder sep subj sep 'm2m_' subj sep subj '.msh'];
elseif (which_pipeline==4)
 roast_mesh_file=dir([subjects_folder sep subj sep '*.msh']);
 [~,idx]=sort([roast_mesh_file.datenum]); %use the oldest msh file assuming that is the original mesh 
 roast_mesh_file=roast_mesh_file(idx);
 msh=mesh_load_gmsh4([subjects_folder sep subj sep roast_mesh_file(1).name]);
 if (length(unique(msh.tetrahedron_regions))>5)
 [tri1, ~, ~]=surftri(mesh_cleaned.node',mesh_cleaned.cell(:,find(mesh_cleaned.field==1))');
 [tri2, ~, ~]=surftri(mesh_cleaned.node',mesh_cleaned.cell(:,find(mesh_cleaned.field==2))');
 [tri3, ~, ~]=surftri(mesh_cleaned.node',mesh_cleaned.cell(:,find(mesh_cleaned.field==3))');
 [tri4, ~, ~]=surftri(mesh_cleaned.node',mesh_cleaned.cell(:,find(mesh_cleaned.field==4))');
 [tri5, ~, ~]=surftri(mesh_cleaned.node',mesh_cleaned.cell(:,find(mesh_cleaned.field==5))');
 tri=[tri1; tri2; tri3; tri4; tri5];
 ind=[ones(1,length(tri1)) 2*ones(1,length(tri2)) 3*ones(1,length(tri3)) 4*ones(1,length(tri4)) 5*ones(1,length(tri5))];
 msh.nodes=mesh_cleaned.node';
 msh.tetrahedra=int32(mesh_cleaned.cell');
 msh.tetrahedron_regions=zeros(length(mesh_cleaned.field),1);
 msh.tetrahedron_regions(:,1)=int32(mesh_cleaned.field);
 msh.triangles=int32(tri);
 msh.triangle_regions=int32(zeros(length(ind),1));
 msh.triangle_regions(:,1)=ind+1000;
 mesh_save_gmsh4(msh, [subjects_folder sep subj sep subj '_reformated' '.msh']);
 tms_opt.fnamehead = [ subjects_folder sep subj sep subj '_reformated' '.msh' ];

 else
  tms_opt.fnamehead = [subjects_folder sep subj sep roast_mesh_file(1).name];
 end
else
 tms_opt.fnamehead = [subjects_folder sep subj sep subj '.msh'];
end
tms_opt.fnamecoil = coil_name;
tms_opt.search_radius = scalp_search_radius;
tms_opt.spatial_resolution = scalp_coil_center_search_grid;
tms_opt.angle_resolution = coil_rotation_discretization;
tms_opt.scalp_normals_smoothing_steps = 50;
tms_opt.open_in_gmsh=0;
if (Optimize_ROI_Efield_Magnitude==0)
 if (TMS_induced_Efield_flows_into_sulcalwall==1)
   target.field=-target.field; %E-Field pointing into sulc wall (equals second part of biphasic pulse)    
 end
end

save([subjects_folder sep subj sep scirunfolder sep 'brain.mat'],'-V6','brain');
save([subjects_folder sep subj sep scirunfolder sep 'target.mat'],'-V6','target');
save([subjects_folder sep subj sep scirunfolder sep 'scalp.mat'],'-V6','scalp');

for i=hairthicknesses
 
 if(i==hairthicknesses(end))
    tms_opt.open_in_gmsh=1; 
 end
   if (Optimize_ROI_Efield_Magnitude==0)
    if (size(target.field,1)==3 && size(target.field,2)==1)
     tms_opt.target_direction = target.field'; 
    else
     tms_opt.target_direction = target.field;   
    end   
  end
 tms_opt.target = target.node';
 if (StartTargetingNavigator && Suppress_Any_Konsole_Input_or_Graphical_Output==0)
 TargetingNavigator('Brain',brain,'Brain_TetCenter',brain_surf_tetcenter,'Scalp',scalp,'ROI_size',tms_opt.target_size,'ROI_center',target.node','ROI_normal',target.field','Scalp_normal',scalp_normal.field');
 input('[TAP] Press enter to continue ', 's');
 if (exist([root_folder sep 'TAP.mat']))
    load([root_folder sep 'TAP.mat'])
    if (isfield(tmp,'target'))
        tms_opt.target=tmp.target;
    end
    if (isfield(tmp,'field'))
        tms_opt.target_direction=tmp.field;
    end
    if (isfield(tmp,'size'))
        tms_opt.target_size=tmp.size;
    end
    if (isfield(tmp,'didt'))
        tms_opt.didt=tmp.didt;
    end 
  end
 end
 already_run=0;
 if (~exist([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i)]))
  tms_opt.pathfem = [subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i)];
  tms_opt.distance=i;
  if (~(size(target.field,1)==1 && size(target.field,2)==3 || size(target.field,1)==3 && size(target.field,2)==1))
      disp('[TAP] Internal error: target direction is not a 3-element vector!');
      break;
  end
  if (~exist([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i)]))
     if (Suppress_Any_Konsole_Input_or_Graphical_Output)
        tms_opt.open_in_gmsh=0;
     end
     run_simnibs(tms_opt);
  else
    disp(['[TAP] For the hairthickness of ' num2str(i) ' mm TAP ran already.']);
    already_run=1;
  end
 end
 if(exist([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i) sep 'target.msh'])~=0)
   msh=mesh_load_gmsh4([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i) '/target.msh']);  
   if (length(msh.element_data)==1)
    roi.node=msh.nodes';
    roi.cell=msh.tetrahedra';
    roi.field=msh.element_data{1}.tetdata;
    roi = clean_mesh(roi.node',roi.cell',roi.field,find(roi.field==0));
    if (exist([subjects_folder sep subj sep scirunfolder sep 'roi.mat' ])==0)
        save([subjects_folder sep subj sep scirunfolder sep 'roi.mat'],'-V6','roi');
    end
 end
 end
 
 if (~already_run)
 log_file=dir([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep '*.log']);
 A=zeros(4,4); %coil orientation in BrainSight = specific transformation matrix
 if (min(size(log_file))==0 && strcmp(coil_name(end-3:end),'.ccd'))
  [optCoil] = optcoil_geo2sci([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep subj '_TMS_optimize_' coil_name(1:end-4) '_coil_pos.geo']);
  save([subjects_folder sep subj sep scirunfolder sep 'optCoilSetup_' target_name  '_HT' sprintf('%.1f',i) '.mat'],'-V6','optCoil');
  %[Arrows, Balls] = coilspace_geo2sci([subj '/' target_name '_HT' sprintf('%.1f',i) '/' 'coil_positions.geo'],'max','none',arrow_lengh);
  %save([subj '/' scirunfolder '/CoilSpace_' target_name '_HT' sprintf('%.1f',i) '.mat'],'-V6','Arrows','Balls');
  canonical_coil.node=canonical_coil.node*1000; %m -> mm
  [R,t] = rigid_transform_3D(canonical_coil.node',optCoil.node');
  A(1:3,1:3)=R;
  A(1:3,4)=t;
 else
  log_file=log_file.name;
  log_file_handle=fopen([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep log_file],'r');
  log_file_text=textscan(log_file_handle,'%s','Delimiter','\n');
  log_file_text=log_file_text{1};
  ind=find(contains(log_file_text,'Best'));
  if (ind+5<=length(log_file_text))
    log_file_text=log_file_text(ind+2:ind+2+3,:);
    log_file_text=strrep(log_file_text,'[','');
    log_file_text=strrep(log_file_text,']','');
    log_file_text=strrep(log_file_text,'  ',' ');
    A(1,:)=sscanf(char(log_file_text(1,:)),'%f %f %f %f');
    A(2,:)=sscanf(char(log_file_text(2,:)),'%f %f %f %f');
    A(3,:)=sscanf(char(log_file_text(3,:)),'%f %f %f %f');
    A(4,:)=sscanf(char(log_file_text(4,:)),'%f %f %f %f');
  end
 end

 A(1:3,2)=-A(1:3,2);
 A(1:3,3)=-A(1:3,3);
 A(4,4)=1;
 
 if (Automatric_Flip_current_direction==1)
   t=A(1:3,4);
   axis_we_measure_angle_from=[0 1 0]';  %y axis pointing anterior-posterior
   axis_we_measure_angle_from=axis_we_measure_angle_from/norm(axis_we_measure_angle_from);
   A(1:3,3)=A(1:3,3)/norm(A(1:3,3));
   tmp=t+axis_we_measure_angle_from;
   proj_point=tmp-dot(tmp-t,A(1:3,3))*A(1:3,3);
   o=proj_point-t;
   o=o / sqrt(sum(o.^2));
   angle=acosd(dot(o, A(1:3,2)));
   if (angle>90 && angle<270)
      A(1:3,2)=-A(1:3,2);  
      Flip_current_direction=1;
   else
      Flip_current_direction=0; 
   end
 end

 n=scalp.node';
 f=scalp.field';
 [mindis index dis] = min_distance( n, coord);
 scalp_normal=f(index,:);
 %check z-axis is correct
 if (dot(A(1:3,3),scalp_normal)<0)
    A(1:3,3)=-A(1:3,3); 
 end
 
 make_brainsight_files(subjects_folder, subj, target_name, A, outputfolder, i, Flip_current_direction, sep, 2);
 
 GM=[];
 if (exist([subjects_folder sep subj sep 'm2m_' subj sep 'gm_fromMesh.nii.gz']))
   GM=load_untouch_nii([ subjects_folder sep subj sep 'm2m_' subj sep 'gm_fromMesh.nii.gz' ]);
   GM=GM.img;
 elseif (exist([subjects_folder sep subj sep 'm2m_' subj sep 'gm.nii.gz']))
   GM=load_untouch_nii([ subjects_folder sep subj sep 'm2m_' subj sep 'gm.nii.gz' ]);
   GM=GM.img;
 elseif (exist([subjects_folder sep subj sep 'm2m_' subj sep subj '_final_contr.nii.gz']))
   GM=load_untouch_nii([ subjects_folder sep subj sep 'm2m_' subj sep subj '_final_contr.nii.gz' ]);
   GM=GM.img;  
 elseif (exist([subjects_folder sep subj sep 'm2m_' subj sep 'final_tissues.nii.gz']))
   GM=load_untouch_nii([subjects_folder sep subj sep 'm2m_' subj sep 'final_tissues.nii.gz']);
   GM=GM.img; 
   GM(find(GM>2))=0;
   GM(find(GM>0 & GM<=2))=1;
 elseif (exist([subjects_folder sep subj sep 'm2m_' subj sep 'surfaces' sep 'cereb_mask.nii.gz']))
   GM=load_untouch_nii( [subjects_folder sep subj sep 'm2m_' subj sep 'surfaces' sep 'cereb_mask.nii.gz'] );
   GM=GM.img;   
   GM(find(GM>0))=1;
 else
   disp(['[TAP] Error: Could not find gray matter segmentation gm.nii.gz or gm_fromMesh.nii.gz. The E-field may be only evaluated in the gray matter and thats done by masking it out in the gray matter. ']);
 end
 
 if (which_pipeline==2) %for mri2mesh TAP can make handknob mask for TMS intensity scaling (Stokes)
  if (exist([root_folder sep 'matlab' sep 'LH_hand_knob.mat'],'file') && ...
      exist([root_folder sep 'matlab' sep 'RH_hand_knob.mat'],'file') && ...
      exist([root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'lh.sphere.reg'],'file') && ...
      exist([root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'rh.sphere.reg'],'file'))
   lh_kn=load([root_folder sep 'matlab' sep 'LH_hand_knob.mat']);
   lh_kn=lh_kn.sphere_lh_hand_knob';
   rh_kn=load([root_folder sep 'matlab' sep 'RH_hand_knob.mat']);
   rh_kn=rh_kn.sphere_rh_hand_knob';
   [sphere_l, ~] = read_surf([ [root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'lh.sphere.reg'] ]);
   [sphere_r, ~] = read_surf([ [root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'rh.sphere.reg'] ]);
   [s_v_l, ind_l]=min_distance(sphere_l, lh_kn);
   [s_v_r, ind_r]=min_distance(sphere_r, rh_kn);

   [pial_l, ~] = read_surf([ [root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'lh.pial'] ]);
   [pial_r, ~] = read_surf([ [root_folder sep subjects_folder sep subj sep 'fs_' subj sep 'surf' sep 'rh.pial'] ]);
   
   scalp_m1_left_distance=min_distance(scalp.node', pial_l(ind_l,:))+hairthicknesses;
   scalp_m1_right_distance=min_distance(scalp.node', pial_r(ind_r,:))+hairthicknesses;
   scalp_target_distance=min_distance(scalp.node', roi.node')+hairthicknesses;
   disp('[TAP] ---- Stokes TMS instensity scaling (assume constant hairthickness) ----');
   disp(['[TAP] adjMT% = ' num2str(Stokes_factor) ' x (coil_target_distance - coil_m1_distance) + ' num2str(MSO)]);
   disp(['[TAP] average/median/min/max scalp distance from ROI= ' num2str(mean(scalp_target_distance)) '/' num2str(median(scalp_target_distance)) '/' num2str(min(scalp_target_distance)) '/' num2str(max(scalp_target_distance)) ' mm']);
   disp(['[TAP] average/median/min/max scalp distance from left Handknob= ' num2str(mean(scalp_m1_left_distance)) '/' num2str(median(scalp_m1_left_distance)) '/' num2str(min(scalp_m1_left_distance)) '/' num2str(max(scalp_m1_left_distance)) ' mm']);
   disp(['[TAP] average/median/min/max scalp distance from right Handknob= ' num2str(mean(scalp_m1_right_distance)) '/' num2str(median(scalp_m1_right_distance)) '/' num2str(min(scalp_m1_right_distance)) '/' num2str(max(scalp_m1_right_distance)) ' mm']);   
   adjustedMSO_left_min= Stokes_factor*(min(scalp_target_distance)-min(scalp_m1_left_distance))+MSO;
   adjustedMSO_left_max= Stokes_factor*(max(scalp_target_distance)-max(scalp_m1_left_distance))+MSO;
   adjustedMSO_left_med= Stokes_factor*(median(scalp_target_distance)-median(scalp_m1_left_distance))+MSO;
   adjustedMSO_left_avr= Stokes_factor*(mean(scalp_target_distance)-mean(scalp_m1_left_distance))+MSO;
    disp(['[TAP] MSO_adjusted for average/median/min/max scalp distance metric= ' num2str(round(adjustedMSO_left_avr)) '/' num2str(round(adjustedMSO_left_med)) '/' num2str(round(adjustedMSO_left_min)) '/' num2str(round(adjustedMSO_left_max)) ' %']);
   disp('[TAP] ---------------------------------------');
   disp('[TAP]')
  end
 end

 if ( (~exist([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep 'res_E.nii.gz'],'file') && Optimize_ROI_Efield_Magnitude==0) || (~exist([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep 'res_normE.nii.gz'],'file') && Optimize_ROI_Efield_Magnitude==1))
  if (which_pipeline==3)
   res=system([simnibs_folder sep 'bin' sep 'msh2nii ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep subj '_TMS_optimize_' coil_name(1:end-4) '.msh '  subjects_folder sep subj sep Original_MRI ' ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files]);
   if (res~=0 && exist(([subjects_folder sep subj sep 'm2m_' subj sep 'T1.nii.gz']),'file'))
     res=system([simnibs_folder sep 'bin' sep 'msh2nii ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep subj '_TMS_optimize_' coil_name(1:end-4) '.msh '  subjects_folder sep subj sep 'm2m_' subj sep 'T1.nii.gz' ' ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files]);
   end
  else
   res=system([simnibs_folder sep 'bin' sep 'msh2nii ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep subj '_TMS_optimize_' coil_name(1:end-4) '.msh '  subjects_folder sep subj sep 'm2m_' subj ' ' subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files]);
  end
 end
 
 if (res==0)
 if (~Use_original_mask_for_Efield_eval)
    Target_nii=load_untouch_nii([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files '_Target.nii.gz']);
 else
     Target_nii=load_untouch_nii(ROI_mask_vol_file);
 end
 
 if (~isempty(GM))
  if (all(size(GM)==size(Target_nii.img)==1))
   gm_ind=intersect(find(Target_nii.img>0),find(GM>0));
   Target_nii.img(:)=0;
   Target_nii.img(gm_ind)=1;
  end
 end
 
 ind=find(Target_nii.img>=Voxelize_Efield_mask_interpolation_thres);
 [x y z] = ind2sub(size(Target_nii.img),find(Target_nii.img>=Voxelize_Efield_mask_interpolation_thres));
 
 if (Optimize_ROI_Efield_Magnitude==1)
  NormE_nii=load_untouch_nii([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files '_normE.nii.gz']);
  NormE_nii=NormE_nii.img;
  ROI_Efield_mag=NormE_nii(ind);
  ROI_Efield_mag=sort(ROI_Efield_mag,'descend');
 
  Emag_value_for_50th_percentile=prctile(ROI_Efield_mag,50);
  Emag_value_for_75th_percentile=prctile(ROI_Efield_mag,75);
  Emag_value_for_90th_percentile=prctile(ROI_Efield_mag,90);
  Emag_value_for_95th_percentile=prctile(ROI_Efield_mag,95);
  Emag_value_for_99th_percentile=prctile(ROI_Efield_mag,99);
 
  if (length(ROI_Efield_mag)>=20)
    Emag_value_for_E20=ROI_Efield_mag(20);
  else
    Emag_value_for_E20=NaN;  
  end
 
  if (length(ROI_Efield_mag)>=50)
    Emag_value_for_E50=ROI_Efield_mag(50);
  else
    Emag_value_for_E50=NaN;  
  end
 
  if (length(ROI_Efield_mag)>=100)
    Emag_value_for_E100=ROI_Efield_mag(100);
  else
    Emag_value_for_E100=NaN;  
  end
  
  if (isnan(Scale_Efield_to_Value) )
    factor=MSO/100*max_didt;
    Emag_value_for_50th_percentile=Emag_value_for_50th_percentile*factor;
    Emag_value_for_75th_percentile=Emag_value_for_75th_percentile*factor;
    Emag_value_for_90th_percentile=Emag_value_for_90th_percentile*factor;
    Emag_value_for_95th_percentile=Emag_value_for_95th_percentile*factor;
    Emag_value_for_99th_percentile=Emag_value_for_99th_percentile*factor;
    Emag_value_for_E100=Emag_value_for_E100*factor;
    Emag_value_for_E50=Emag_value_for_E50*factor;
    Emag_value_for_E20=Emag_value_for_E20*factor;
    outputformat=['%.' num2str(Efield_display_decimal_places) 'f'];
    disp(['[TAP] ROI E-field magnitude (50th percentile) = ' sprintf(outputformat, Emag_value_for_50th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (75th percentile) = ' sprintf(outputformat, Emag_value_for_75th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (90th percentile) = ' sprintf(outputformat, Emag_value_for_90th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (95th percentile) = ' sprintf(outputformat, Emag_value_for_95th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (99th percentile) = ' sprintf(outputformat, Emag_value_for_99th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (E100) = ' sprintf(outputformat, Emag_value_for_E100) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (E50) = ' sprintf(outputformat, Emag_value_for_E50) ' V/m']);
    disp(['[TAP] ROI E-field magnitude (E20) = ' sprintf(outputformat, Emag_value_for_E20) ' V/m']);
  elseif isnumeric(Scale_Efield_to_Value)
    Scaling_value_for_50th_percentile=Scale_Efield_to_Value/Emag_value_for_50th_percentile;
    Scaling_value_for_75th_percentile=Scale_Efield_to_Value/Emag_value_for_75th_percentile;
    Scaling_value_for_90th_percentile=Scale_Efield_to_Value/Emag_value_for_90th_percentile;
    Scaling_value_for_95th_percentile=Scale_Efield_to_Value/Emag_value_for_95th_percentile;
    Scaling_value_for_99th_percentile=Scale_Efield_to_Value/Emag_value_for_99th_percentile;
    Scaling_value_for_E100=Scale_Efield_to_Value/Emag_value_for_E100;
    Scaling_value_for_E50=Scale_Efield_to_Value/Emag_value_for_E50;
    Scaling_value_for_E20=Scale_Efield_to_Value/Emag_value_for_E20; 
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 50th percentile of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_50th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_50th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 75th percentile of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_75th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_75th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 90th percentile of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_90th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_90th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 95th percentile of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_95th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_95th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 99th percentile of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_99th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_99th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E100 of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_E100/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E100)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E50 of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_E50/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E50)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E20 of ROI E-field magnitude is %MSO = ' num2str(round(Scaling_value_for_E20/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E20)) ' A/' sprintf('\x3BC') 's)']);
  end

 else
  vecE_nii=load_untouch_nii([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep prefix_for_voxelized_nii_files '_E.nii.gz']);
  vecE_nii=vecE_nii.img;
 
  ROI_Efield_vec=zeros(length(x),1);
  if (length(x)==length(y) && length(y)==length(z))
   for j=1:length(x)
    ROI_Efield_vec(j,1)=dot(target.field,squeeze(vecE_nii(x(j),y(j),z(j),:)));
   end
  end
  
  ROI_Efield_vec=sort(ROI_Efield_vec,'descend');
  
  if (length(ROI_Efield_vec)==1)
   Evec_value_for_50th_percentile=NaN;
   Evec_value_for_75th_percentile=NaN;
   Evec_value_for_90th_percentile=NaN;
   Evec_value_for_95th_percentile=NaN;
   Evec_value_for_99th_percentile=NaN;
  else
   Evec_value_for_50th_percentile=prctile(ROI_Efield_vec,50);
   Evec_value_for_75th_percentile=prctile(ROI_Efield_vec,75);
   Evec_value_for_90th_percentile=prctile(ROI_Efield_vec,90);
   Evec_value_for_95th_percentile=prctile(ROI_Efield_vec,95);
   Evec_value_for_99th_percentile=prctile(ROI_Efield_vec,99);
  end
  
  if (length(ROI_Efield_vec)>=20)
    Evec_value_for_E20=ROI_Efield_vec(20);
  else
    Evec_value_for_E20=NaN;  
  end
 
  if (length(ROI_Efield_vec)>=50)
    Evec_value_for_E50=ROI_Efield_vec(50);
  else
    Evec_value_for_E50=NaN;  
  end
 
  if (length(ROI_Efield_vec)>=100)
    Evec_value_for_E100=ROI_Efield_vec(100);
  else
    Evec_value_for_E100=NaN;  
  end
  
  if (isnan(Scale_Efield_to_Value) && max_didt~=1 && max_didt~=0 && ~isnan(max_didt))
    factor=MSO/100*max_didt;
    Evec_value_for_50th_percentile=Evec_value_for_50th_percentile*factor;
    Evec_value_for_75th_percentile=Evec_value_for_75th_percentile*factor;
    Evec_value_for_90th_percentile=Evec_value_for_90th_percentile*factor;
    Evec_value_for_95th_percentile=Evec_value_for_95th_percentile*factor;
    Evec_value_for_99th_percentile=Evec_value_for_99th_percentile*factor;
    Evec_value_for_E100=Evec_value_for_E100*factor;
    Evec_value_for_E50=Evec_value_for_E50*factor;
    Evec_value_for_E20=Evec_value_for_E20*factor;
    outputformat=['%.' num2str(Efield_display_decimal_places) 'f'];
    disp(['[TAP] ROI E-field normal (50th percentile) = ' sprintf(outputformat, Evec_value_for_50th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field normal (75th percentile) = ' sprintf(outputformat, Evec_value_for_75th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field normal (90th percentile) = ' sprintf(outputformat, Evec_value_for_90th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field normal (95th percentile) = ' sprintf(outputformat, Evec_value_for_95th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field normal (99th percentile) = ' sprintf(outputformat, Evec_value_for_99th_percentile) ' V/m']);
    disp(['[TAP] ROI E-field normal (E100) = ' sprintf(outputformat, Evec_value_for_E100) ' V/m']);
    disp(['[TAP] ROI E-field normal (E50) = ' sprintf(outputformat, Evec_value_for_E50) ' V/m']);
    disp(['[TAP] ROI E-field normal (E20) = ' sprintf(outputformat, Evec_value_for_E20) ' V/m']);
  elseif (isnumeric(Scale_Efield_to_Value) && max_didt~=1 && max_didt~=0 && ~isnan(max_didt))
    Scaling_value_for_50th_percentile=Scale_Efield_to_Value/Evec_value_for_50th_percentile;
    Scaling_value_for_75th_percentile=Scale_Efield_to_Value/Evec_value_for_75th_percentile;
    Scaling_value_for_90th_percentile=Scale_Efield_to_Value/Evec_value_for_90th_percentile;
    Scaling_value_for_95th_percentile=Scale_Efield_to_Value/Evec_value_for_95th_percentile;
    Scaling_value_for_99th_percentile=Scale_Efield_to_Value/Evec_value_for_99th_percentile;
    Scaling_value_for_E100=Scale_Efield_to_Value/Evec_value_for_E100;
    Scaling_value_for_E50=Scale_Efield_to_Value/Evec_value_for_E50;
    Scaling_value_for_E20=Scale_Efield_to_Value/Evec_value_for_E20; 
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 50th percentile of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_50th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_50th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 75th percentile of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_75th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_75th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 90th percentile of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_90th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_90th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 95th percentile of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_95th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_95th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for 99th percentile of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_99th_percentile/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_99th_percentile)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E100 of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_E100/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E100)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E50 of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_E50/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E50)) ' A/' sprintf('\x3BC') 's)']);
    disp(['[TAP] To reach ' num2str(Scale_Efield_to_Value) ' V/m, the needed scaling for E20 of ROI E-field normal is %MSO = ' num2str(round(Scaling_value_for_E20/max_didt*100)) ' (didt = ' num2str(round(Scaling_value_for_E20)) ' A/' sprintf('\x3BC') 's)']);
  end
 end
 end
 system(char(['rm -rf ' subjects_folder sep subj sep' target_name '_HT' sprintf('%1.1f',i)]));
 end
 disp(['[TAP] --- computations for hairthickness ' num2str(i) ' mm completed ---']);
end

disp(['[TAP] All computations done for subject ' subj '.']);