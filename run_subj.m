clear all;

Flip_current_direction=1; %1=true, 0=no
Angle_tolerance_flipcoilzaxis=30;

subj='24372';
target_name='PMC';
mask='24372_PMC.nii.gz';
outputfolder_deliver_to_norren='PMC_4Noreen';
scirunfolder='PMC_scirun';

simnibs_folder='/Users/rothermichk17/Applications/SimNIBS-3.2/simnibs/';
root_folder='/Users/rothermichk17/Pipeline';
freesurfer_matlab_folder='/Applications/freesurfer/7.1.0//matlab/';

normal_search_radius=4;
hairthicknesses=0.5:0.5:3.0;%:0.5:7.5; %in mm
scalp_search_radius=25; %mm radius around scalp-prpjected point
scalp_coil_center_search_grid=1; %1mm Manhattan distance
coil_rotation_discretization = 1; %1 degree
arrow_lengh=0.4;
coil_name='Medtronic_MCF_B65.ccd';
addpath([simnibs_folder '/matlab']);
addpath(char([root_folder '/matlab']));
addpath(char(freesurfer_matlab_folder));

coord = get_target(subj,mask);
[brain, scalp, target] = get_target_normal(subj, coord, [1 2], normal_search_radius);
brainsight_brain_location = convert_coord_simnibs_2_brainsight(subj, mask, target.node');
disp([ target_name ' BRAINSIGHT COORDINATE: ' sprintf('%.3f %.3f %.3f',brainsight_brain_location(1),brainsight_brain_location(2),brainsight_brain_location(3))]);
if (~exist([subj '/' scirunfolder])) 
  mkdir([subj '/' scirunfolder]);
end
save([subj '/' scirunfolder '/brain.mat'],'-V6','brain');
save([subj '/' scirunfolder '/target.mat'],'-V6','target');
save([subj '/' scirunfolder '/scalp.mat'],'-V6','scalp');

tms_opt = opt_struct('TMSoptimize');
tms_opt.fnamehead = [subj '/' subj '_bin.msh'];
tms_opt.fnamecoil = coil_name;
tms_opt.method = 'ADM';
tms_opt.search_radius = scalp_search_radius;
tms_opt.spatial_resolution = scalp_coil_center_search_grid;
tms_opt.search_angle = 360;
tms_opt.angle_resolution = coil_rotation_discretization;
tms_opt.scalp_normals_smoothing_steps = 123;
tms_opt.open_in_gmsh=0;
target.field=-target.field; %E-Field pointing into sulc wall (second part of biphasic pulse)
for i=hairthicknesses
 if(i==hairthicknesses(end))
    tms_opt.open_in_gmsh=1; 
 end
 if (~exist([subj '/' target_name '_HT' sprintf('%.1f',i)]))
  tms_opt.pathfem = [subj '/' target_name '_HT' sprintf('%.1f',i)];
  tms_opt.distance=i;
  tms_opt.target = target.node';
  if (~(size(target.field,1)==1 && size(target.field,2)==3 || size(target.field,1)==3 && size(target.field,2)==1))
      disp('Internal error: target direction is not a 3-element vector!');
      break;
  end  
  if (size(target.field,1)==3 && size(target.field,2)==1)
    tms_opt.target_direction = -target.field'; 
  else
    tms_opt.target_direction = -target.field;   
  end
  run_simnibs(tms_opt);
 end
 if(exist([subj '/' target_name '_HT' sprintf('%.1f',i) '/target.msh'])~=0)
   msh=mesh_load_gmsh4([subj '/' target_name '_HT' sprintf('%.1f',i) '/target.msh']);  
   if (length(msh.element_data)==2)
    roi.node=msh.nodes';
    roi.cell=msh.tetrahedra';
    roi.field=msh.element_data{2}.tetdata;
    roi = clean_mesh(roi.node',roi.cell',roi.field,find(roi.field==0));
    if (exist([subj '/' scirunfolder '/roi.mat' ])==0)
        save([subj '/' scirunfolder '/roi.mat'],'-V6','roi');
    end
 end
 end
 
[optCoil] = optcoil_geo2sci([subj '/' target_name  '_HT' sprintf('%.1f',i) '/' subj '_bin' '_TMS_optimize_' coil_name(1:end-4) '_coil_pos.geo']);
save([subj '/' scirunfolder '/optCoilSetup_' target_name  '_HT' sprintf('%.1f',i) '.mat'],'-V6','optCoil');
%[Arrows, Balls] = coilspace_geo2sci([subj '/' target_name '_HT' sprintf('%.1f',i) '/' 'coil_positions.geo'],'max','none',arrow_lengh);
%save([subj '/' scirunfolder '/CoilSpace_' target_name '_HT' sprintf('%.1f',i) '.mat'],'-V6','Arrows','Balls');
canonical_coil=load_coil([simnibs_folder '/simnibs_env/lib/python3.7/site-packages/simnibs/ccd-files/' coil_name]);%load coil definition so a transfermatrix to the optimal setup can be computed
canonical_coil.node=canonical_coil.node*1000; %m -> mm
[R,t] = rigid_transform_3D(canonical_coil.node',optCoil.node');
brainsight = convert_coord_simnibs_2_brainsight(subj, mask, t');
A=zeros(4,4); %coil orientation in BrainSight = specific transformation matrix - empirically determined and visually verified
R(:,3)=-R(:,3);
A(1:3,1)=cross(R(:,2),R(:,3));
A(1:3,2)=R(:,2);
A(1:3,3)=R(:,3);
A(3,1:2)=-A(3,1:2); 
A(1:3,4)=brainsight;
A(4,4)=1;

if (Flip_current_direction==1)
    A(1:3,2)=-A(1:3,2); 
end

%check if coil's z-axis need to be flipped  
[mindis index dis] = min_distance(scalp.node', t');
scalp_normal=scalp.field(:,index)';
scalp_normal=scalp_normal/norm(scalp_normal);
coil_z_axis=A(1:3,3)';
coil_z_axis=coil_z_axis/norm(coil_z_axis);
disp(['angle: ' num2str(acosd(dot(scalp_normal, coil_z_axis)))]);
 if (abs(acosd(dot(scalp_normal, coil_z_axis))-180) < Angle_tolerance_flipcoilzaxis)  %(Flip_coil_zaxis==1)
    A(1:3,3)=-A(1:3,3);  
 end

make_brainsight_files(subj, target_name, A, outputfolder_deliver_to_norren, i, Flip_current_direction);
disp(['--- computations for hairthickness ' num2str(i) ' mm completed ---']);

end

disp(['All computations for ' subj ' completed.']);
 
 

