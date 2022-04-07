clear all;

TMS_COIL_ANGLE_TO_MIDLINE = 45; %degree
HairThinkness = 10; %mm

if (ispc)
    sep='\';
elseif (ismac || isunix)
    sep='/';
end

subj='Ernie_headreco';
subjects_folder='subjects';
target_name='DLPFC';
mask='DLPFC_mask_TMS_target.nii.gz';%'M1_mask_TMS_target.nii.gz';%'fiducial_AIL.nii.gz';%'fiducial_AIL.nii.gz';%'gm_fromMesh.nii.gz';%'fiducial_AIL.nii.gz';%M1_mask_TMS_target.nii.gz';%'fiducials_AIL.nii.gz';%'fiducials.nii.gz';%'gm_fromMesh.nii.gz';%'fiducials_AIL.nii.gz';%'fiducials.nii.gz';%'ernie_masks_contr.nii.gz';%'fiducial_points.nii.gz';%'gm_fromMesh.nii.gz';%'M1_mask_opt.nii.gz';%'gm_fromMesh.nii.gz';%'ernie_masks_contr.nii.gz';%'gm_fromMesh.nii.gz';%'ernie_masks_contr.nii.gz';%'fiducial_points.nii.gz';%'M1_mask_opt.nii.gz';%exp/24372_PMC_writingactive_firstlevel_zstat1_T1w.nii.gz';%'gm_only.nii.gz';%'exp/PMC_RUA.nii.gz';%'exp/24372_PMC_LDB.nii.gz'; %'24372_PMC.nii.gz'; %'gm_only.nii.gz';
outputfolder='DLPFC_TMSsim';
scirunfolder='DLPFC_scirun';
SimNIBS_MRI=[subj '_T1fs_conform.nii.gz'];

simnibs_folder='/Users/rothermichk17/Applications/SimNIBS-3.2/';
root_folder='/Users/rothermichk17/Pipeline_git';
freesurfer_matlab_folder='/Applications/freesurfer/7.1.0//matlab/';

coil_name='Medtronic_MCF_B65.ccd';
addpath(char([simnibs_folder sep 'matlab']));
addpath(char([root_folder sep 'matlab']));
addpath(char(freesurfer_matlab_folder));

[brain, scalp, head_model_mesh, brain_vert_neigh]=get_meshSurfs(subjects_folder, subj, [1 2], sep);
[coord, which_pipeline] = get_target(subjects_folder,subj,mask,SimNIBS_MRI, sep, brain);

scalp_nodes=scalp.node';
[~, index, ~] = min_distance(scalp_nodes, coord);
spp=scalp_nodes(index,:);
zdir=scalp.field(:,index)';

v=[0 1 0]; %ant-post direction
ydir=cross(v, zdir);
ydir=ydir/norm(ydir);
xdir=cross(ydir, zdir);
xdir=xdir/norm(xdir);

u=zdir;
ux=u(1);
uy=u(2);
uz=u(3);
c=cosd(TMS_COIL_ANGLE_TO_MIDLINE);
s=sind(TMS_COIL_ANGLE_TO_MIDLINE);
R=[c+ux^2*(1-c) ux*uy*(1-c)-uz*s ux*uz*(1-c)+uy*s; 
     uy*ux*(1-c)+uz*s c+uy^2*(1-c) uy*uz*(1-c)-ux*s; 
     uz*ux*(1-c)-uy*s uz*uy*(1-c)+ux*s c+uz^2*(1-c)
     ]; 

xdir=R*xdir';
ydir=R*ydir';

%for visualization in scirun
Field1.node=spp';
Field1.field=xdir;
Field2.node=spp';
Field2.field=ydir;
Field3.node=spp';
Field3.field=zdir';

T=zeros(4,4);
T(1:3,1)=xdir;
T(1:3,2)=ydir;
T(1:3,3)=-zdir;
tmp=spp+HairThinkness.*zdir;
T(1:3,4)=tmp';
T(4,4)=1;

Flip_current_direction=0;

S = sim_struct('SESSION');
S.fnamehead = [subjects_folder sep subj sep 'ernie_headreco_no-conform.msh']; % head mesh
S.pathfem = [ subjects_folder sep subj sep outputfolder];
S.poslist{1} = sim_struct('TMSLIST');
S.poslist{1}.fnamecoil = 'Magstim_70mm_Fig8.ccd'; % Choose a coil in the ccd-files folder
S.poslist{1}.anisotropy_type='scalar';
S.poslist{1}.aniso_maxcond=2;
S.poslist{1}.aniso_maxratio=10;
S.poslist{1}.pos(1).distance = 0; % 4 mm distance from coil surface to head surface
S.poslist{1}.pos(1).didt = 1e6; % 4 mm distance from coil surface to head surface
S.poslist{1}.pos(1).matsimnibs=T;
run_simnibs(S);

A=T;
make_brainsight_files(subjects_folder, subj, target_name, A, outputfolder, HairThinkness, Flip_current_direction, sep, 2);
 
