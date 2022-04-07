clear all;

if (ispc)
    sep='\';
elseif (ismac || isunix)
    sep='/';
end

Automatric_Flip_current_direction= 1;

subj='Ernie_headreco';
subjects_folder='subjects';
target_name='M1';
mask='M1_mask_TMS_target.nii.gz';
scirunfolder='M1_scirun';
SimNIBS_MRI=[subj '_T1fs_conform.nii.gz'];
tms_opt.target_size=3; %radius of ROI 

simnibs_folder='/Users/rothermichk17/Applications/SimNIBS-3.2/';
root_folder='/Users/rothermichk17/Pipeline_git';
freesurfer_matlab_folder='/Applications/freesurfer/7.1.0//matlab/';

hairthicknesses=0.5;%:0.5:3.0;%:0.5:7.5; %in mm
scalp_search_radius=25; %mm radius around scalp-prpjected point
scalp_coil_center_search_grid=1; %1mm Manhattan distance
coil_rotation_discretization = 1; %1 degree
arrow_lengh=0.4;
coil_name='Medtronic_MCF_B65.ccd';
addpath(char([simnibs_folder sep 'matlab']));
addpath(char([root_folder sep 'matlab']));
addpath(char(freesurfer_matlab_folder));

[brain, scalp, head_model_mesh, brain_vert_neigh]=get_meshSurfs(subjects_folder, subj, [1 2], sep);
[coord, which_pipeline] = get_target(subjects_folder,subj,mask,SimNIBS_MRI, sep, brain);
[target] = get_target_normal(subjects_folder, sep, subj, coord, brain, scalp, head_model_mesh, which_pipeline, brain_vert_neigh);

if (~exist([subjects_folder sep subj sep scirunfolder])) 
  mkdir([subjects_folder sep subj sep scirunfolder]);
end
save([subjects_folder sep subj sep scirunfolder '/brain.mat'],'-V6','brain');
save([subjects_folder sep subj sep scirunfolder '/target.mat'],'-V6','target');
save([subjects_folder sep subj sep scirunfolder '/scalp.mat'],'-V6','scalp');

if (~exist([subj '/' outputfolder ])) 
  mkdir([subjects_folder sep subj sep outputfolder ]);
end

if (which_pipeline==1)
  stlwrite([subjects_folder sep subj sep outputfolder sep subj '_LoadAsNIFTI_brain.stl'],brain.face',brain.node');
  stlwrite([subjects_folder sep subj sep outputfolder sep subj '_LoadAsNIFTI_scalp.stl'],scalp.face',scalp.node'); 
end

tms_opt = opt_struct('TMSoptimize');
tms_opt.fnamehead = [subjects_folder sep subj sep subj '.msh'];
tms_opt.fnamecoil = coil_name;
tms_opt.method = 'ADM';
tms_opt.search_radius = scalp_search_radius;
tms_opt.spatial_resolution = scalp_coil_center_search_grid;
tms_opt.search_angle = 360;
tms_opt.angle_resolution = coil_rotation_discretization;
tms_opt.scalp_normals_smoothing_steps = 50;
tms_opt.open_in_gmsh=0;
tms_opt.target_size=3;
target.field=-target.field; %E-Field pointing into sulc wall (second part of biphasic pulse)

for i=hairthicknesses
 if(i==hairthicknesses(end))
    tms_opt.open_in_gmsh=1; 
 end
 if (~exist([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i)]))
  tms_opt.pathfem = [subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i)];
  tms_opt.distance=i;
  tms_opt.target = target.node';
  if (~(size(target.field,1)==1 && size(target.field,2)==3 || size(target.field,1)==3 && size(target.field,2)==1))
      disp('Internal error: target direction is not a 3-element vector!');
      break;
  end  
  if (size(target.field,1)==3 && size(target.field,2)==1)
    tms_opt.target_direction = target.field'; 
  else
    tms_opt.target_direction = target.field;   
  end
  run_simnibs(tms_opt);
 end
 if(exist([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i) sep 'target.msh'])~=0)
   msh=mesh_load_gmsh4([subjects_folder sep subj sep target_name '_HT' sprintf('%.1f',i) '/target.msh']);  
   if (length(msh.element_data)==2)
    roi.node=msh.nodes';
    roi.cell=msh.tetrahedra';
    roi.field=msh.element_data{2}.tetdata;
    roi = clean_mesh(roi.node',roi.cell',roi.field,find(roi.field==0));
    if (exist([subjects_folder sep subj sep scirunfolder sep 'roi.mat' ])==0)
        save([subjects_folder sep subj sep scirunfolder sep 'roi.mat'],'-V6','roi');
    end
 end
 end
 
 [optCoil] = optcoil_geo2sci([subjects_folder sep subj sep target_name  '_HT' sprintf('%.1f',i) sep subj '_TMS_optimize_' coil_name(1:end-4) '_coil_pos.geo']);
 save([subjects_folder sep subj sep scirunfolder sep 'optCoilSetup_' target_name  '_HT' sprintf('%.1f',i) '.mat'],'-V6','optCoil');
 %[Arrows, Balls] = coilspace_geo2sci([subj '/' target_name '_HT' sprintf('%.1f',i) '/' 'coil_positions.geo'],'max','none',arrow_lengh);
 %save([subj '/' scirunfolder '/CoilSpace_' target_name '_HT' sprintf('%.1f',i) '.mat'],'-V6','Arrows','Balls');
 canonical_coil=load_coil([simnibs_folder '/simnibs_env/lib/python3.7/site-packages/simnibs/ccd-files/' coil_name]);%load coil definition so a transfermatrix to the optimal setup can be computed
 canonical_coil.node=canonical_coil.node*1000; %m -> mm
 [R,t] = rigid_transform_3D(canonical_coil.node',optCoil.node');
 
 A=zeros(4,4); %coil orientation in BrainSight = specific transformation matrix - empirically determined and visually verified
 A(1:3,1:3)=R;
 A(1:3,2)=-A(1:3,2);
 A(1:3,3)=-A(1:3,3);
 A(1:3,4)=t;
 A(4,4)=1;
 
 if (Automatric_Flip_current_direction==1)
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

 make_brainsight_files(subjects_folder, subj, target_name, A, outputfolder, i, Flip_current_direction, sep, 2);
 system(char(['rm -rf ' subjects_folder subj sep' target_name '_HT' sprintf('%1.1f',i)]));
 disp(['--- computations for hairthickness ' num2str(i) ' mm completed ---']);
end

disp(['All computations for ' subj ' completed.']);



