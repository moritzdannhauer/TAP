function [coord, ROI_size, which_pipeline, ROI_mask_vol_file, simnibsModel] = get_target(subjects_folder, subject, mask, mri, sep, brain, Original_MRI, fsl_path)
ROI_size=0;
coord=NaN;
epsi=1e-8;
which_pipeline=0;
alignment_voxel_stepsize=1000; %use every 1000th voxel of gm mask to match brain surface
fine_alignment_voxel_stepsize=500;
fine_alignment_number_voxel_to_match=150;
tri_neighborhood=30; % [mm]
voxelspace_to_mesh_moving_range=-2:1:2;
alignment_range=max(abs(voxelspace_to_mesh_moving_range))+1;
charm_gm_mask='cereb_mask.nii.gz';
gm_mask='gm.nii.gz';
roast_masks='*masks.nii.gz';
resampled_mask='ResampledMask.nii.gz';
tap_extracted_org_brain_file_name='TAP_extracted_org_brain.nii.gz';
tap_extracted_simnibs_brain_file_name='TAP_extracted_simnibs_brain.nii.gz';
tap_tmp_file='TAP_tmp.txt';
tap_tmp_nii_file='TAP_tmp.nii.gz';
tap_trans_nii_file='TAP_affine_transformation_matrix.txt';

if (nargin~=8 || ~exist('subjects_folder','var') || ~exist('subject','var') || ~exist('mask','var') || ~exist('mri','var') || ~exist('sep','var') || ~exist('brain','var') || ~exist('brain','var') || ~exist('Original_MRI','var') || ~exist('fsl_path','var'))
  disp('[TAP] TAP''s function get_target was not called with the right number of parameter. ');
  return ;
end

if exist(([subjects_folder sep subject sep 'm2m_' subject sep 'headreco_log.html']),'file')
  disp(['[TAP] ' subject ' processed with headreco pipeline. ' ]);
  which_pipeline=1;
elseif exist(([subjects_folder sep subject sep 'm2m_' subject sep 'mri2mesh_log.html']),'file')
   disp(['[TAP] ' subject ' processed with mri2mesh pipeline. ' ]);   
   which_pipeline=2;
elseif  exist(([subjects_folder sep subject sep 'm2m_' subject sep 'charm_log.html']),'file')
   disp(['[TAP] ' subject ' processed with charm pipeline. ' ]);   
   which_pipeline=3; 
else 
%must be ROAST then
   disp(['[TAP] No information was found for ' subject ' so TAP assumes it is from ROAST. ' ]);   
   which_pipeline=4; 
end

ROI_mask_vol_file=[subjects_folder sep subject sep mask];
if ~exist((ROI_mask_vol_file),'file')
  error(['[TAP] ERROR: File '  [subject sep mask] ' does not exist. ' ]);
end

[TMaskMRI, p, Mask_MRI_CS, dimMaskMRI]=load_nii_file(ROI_mask_vol_file, eps);

if ~exist(([subjects_folder sep subject sep mri]),'file')
  if exist(([subjects_folder sep subject sep 'm2m_' subject sep 'T1.nii.gz']),'file')
      mri=['m2m_' subject sep 'T1.nii.gz'];
  else
    error(['[TAP] ERROR: File '  [subject sep mri] ' does not exist. ' ]);
  end
end

[TsimNIBSMRI, ~, SimNIBS_MRI_CS, dimSimnibsMRI]=load_nii_file([subjects_folder sep subject sep mri],0);

if (~strcmp(Mask_MRI_CS,SimNIBS_MRI_CS) || any(dimSimnibsMRI~=dimMaskMRI))
   disp(['[TAP] Ok, the mask and SimNIBS-generated MRI you provided are in different spaces. TAP will assume the mask was registered to the original MRI and try to register it to the SimNIBS-generated MRI.']); 
    if (exist(([subjects_folder sep subject sep Original_MRI]),'file'))
     [TOrgMRI, ~, Org_MRI_CS, dimOrgMRI]=load_nii_file([subjects_folder sep subject sep Original_MRI],1-eps);
%      if (~strcmp(Org_MRI_CS,Mask_MRI_CS) || any(dimOrgMRI~=dimMaskMRI))
%         disp(['[TAP] Also the original MRI and the mask are in different space, you need to register the mask to the original MRI to run TAP, like:']);
%         disp(['[TAP] ' fsl_path sep 'bin' sep 'flirt -in ' ROI_mask_vol_file ' -ref ' subjects_folder sep subject sep Original_MRI ' -out ' resampled_mask ' -applyxfm']);        
%         disp(['[TAP] ' fsl_path sep 'bin' sep 'flirt -in skullstrippedMRI_the_mask_is_from.nii.gz -ref skullstrippedSimNIBSMRI_T1fs_conform.nii.gz -out RegisteredMRI.nii.gz -omat AffineTransformationMatrix.txt']);
%         disp(['[TAP] ' fsl_path sep 'bin' sep 'flirt -in ' subjects_folder sep subject sep mask  '  -ref RegisteredMRI.nii.gz -applyxfm -init AffineTransformationMatrix.txt -out ' subjects_folder sep subject sep 'trans_' mask]);
%         disp(['[TAP] And then use ' subjects_folder sep subject sep 'trans_' mask ' as input for TAP']);
%         return;
%      end
     [fsl1, ~] =system([fsl_path sep 'bet' ' > ' subjects_folder sep subject sep tap_tmp_file]);  %check if FSL's bet is installed
     [fsl2, ~] =system([fsl_path sep 'flirt' ' > ' subjects_folder sep subject sep tap_tmp_file]); 
     if (~fsl1 || ~fsl2 || ~exist([subjects_folder sep subject sep tap_tmp_file],'file'))
       disp('[TAP] FSL bet is not accessable from MATLAB. Lets try without FSL then.');
     else
       disp('[TAP] Running FSL/flirt - this may take a few minutes ');  
       setenv('FSLDIR',char(fsl_path));
       setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); 
       [org_brain_extracting, ~]=system(char([fsl_path sep 'bin' sep 'bet ' subjects_folder sep subject sep Original_MRI ' ' subjects_folder sep subject sep tap_extracted_org_brain_file_name ' > ' subjects_folder sep subject sep tap_tmp_file]));
       [simnibs_brain_extracting, ~]=system(char([fsl_path sep 'bin' sep 'bet ' subjects_folder sep subject sep mri  ' ' subjects_folder sep subject sep tap_extracted_simnibs_brain_file_name ' > ' subjects_folder sep subject sep tap_tmp_file]));
       if (org_brain_extracting==0 && simnibs_brain_extracting==0)
          [register_brains, ~] = system([fsl_path sep 'bin' sep 'flirt -in ' subjects_folder sep subject sep tap_extracted_org_brain_file_name  '  -ref ' subjects_folder sep subject sep tap_extracted_simnibs_brain_file_name ' -out ' subjects_folder sep subject sep tap_tmp_nii_file ' -omat ' subjects_folder sep subject sep tap_trans_nii_file]);
          if (any(dimOrgMRI~=dimMaskMRI) || ~strcmp(Org_MRI_CS,Mask_MRI_CS))
              disp(['[TAP] Also the original MRI and the mask are in different space, TAP is going to try to correct that, this can FAIL, so, ']);
              disp(['[TAP] the FSL-Flirt call is executed but PLEASE overlay the result ' mask ' and ' Original_MRI ' in your MRI viewer to make sure - there no way for TAP to check automatically!!']);    
              [resample_mask_to_org_mri, ~] = system([fsl_path sep 'bin' sep 'flirt -in ' subjects_folder sep subject sep mask  '  -ref ' subjects_folder sep subject sep Original_MRI ' -out ' subjects_folder sep subject sep 'trans_' mask ' -applyxfm']);
              if (resample_mask_to_org_mri==0)
                 mask=['trans_' mask];
                 if (isempty(strfind(mask,'.nii.gz')))
                    if (~isempty(strfind(mask,'.nii'))) 
                        mask=[mask '.gz'];
                    else
                        mask=[mask '.nii.gz'];
                    end
                 end
                 disp(['[TAP] The FSL-Flirt call was executed but please overlay ' mask ' and ' Original_MRI ' in your MRI viewer to make sure - there no way for TAP to check!!']);
              else
                  disp(['[TAP] TAP tried its best but the FSL-Flirt command errored out. Aborting this TAP run!']);
                  return;
              end
          end
          [move_target_to_simnibsbrains, ~] = system([fsl_path sep 'bin' sep 'flirt -in ' subjects_folder sep subject sep mask  '  -ref ' subjects_folder sep subject sep tap_tmp_nii_file ' -applyxfm -init ' subjects_folder sep subject sep tap_trans_nii_file ' -out ' subjects_folder sep subject sep 'trans_' mask]);
       end
       if (exist([subjects_folder sep subject sep 'trans_' mask],'file'))
        str=[subjects_folder sep subject sep 'trans_' mask];
        nr_voxels=size(p,1);
        [TMaskMRI, p, Mask_MRI_CS, ~, Voxels]=load_nii_file(str,eps);
        tmp_ind = sub2ind(size(Voxels),p(:,1),p(:,2),p(:,3));
        Voxels=Voxels(tmp_ind);
        [~,tmp_ind] = sort(Voxels,'descend');
        p=p(tmp_ind(1:nr_voxels),:);
        ROI_mask_vol_file=str;
       else
        disp(['[TAB] cound not find registered mask: ' str]);
        coord=NaN;
        return ;
       end
     end
     if (ispc)
       system(['del /f ' subjects_folder sep subject sep tap_tmp_file]);
       system(['del /f ' subjects_folder sep subject sep tap_tmp_nii_file]);
       system(['del /f ' subjects_folder sep subject sep tap_trans_nii_file]);
       system(['del /f ' subjects_folder sep subject sep tap_extracted_org_brain_file_name]);
       system(['del /f ' subjects_folder sep subject sep tap_extracted_simnibs_brain_file_name]);
     elseif (ismac || isunix)
       system(['rm -f ' subjects_folder sep subject sep tap_tmp_file]);
       system(['rm -f ' subjects_folder sep subject sep tap_tmp_nii_file]);
       system(['rm -f ' subjects_folder sep subject sep tap_trans_nii_file]);
       system(['rm -f ' subjects_folder sep subject sep tap_extracted_org_brain_file_name]);
       system(['rm -f ' subjects_folder sep subject sep tap_extracted_simnibs_brain_file_name]);
     end
    else
     disp(['[TAP] TAP cannot find the original MRI: ' Original_MRI '. Check the path to the original MRI.']);
     coord=NaN;
     return ;
    end
end

%fine alignment MRI and mesh
tmp1=exist([subjects_folder sep subject sep 'm2m_' subject sep gm_mask],'file');
tmp2=0;
if (which_pipeline==3)
  tmp2=exist([subjects_folder sep subject sep 'm2m_' subject sep 'surfaces' sep charm_gm_mask],'file');
end
if (which_pipeline==4)
  roast_mask_file=dir([subjects_folder sep subject sep roast_masks]);
  [~,idx]=sort([roast_mask_file.datenum]); %use the oldest msh file assuming that is the original mesh 
  roast_mask_file=roast_mask_file(idx);
  tmp2=exist([subjects_folder sep subject sep roast_mask_file(1).name],'file');
end

if (tmp1~=0 || tmp2~=0) 
    if (tmp1~=0)
      disp(['[TAP] Found ' gm_mask ' and now refine alignment (the following steps may take a few minutes, have patience) ']);
      [~, ~, ~, ~, GM_fromMesh] = load_nii_file([subjects_folder sep subject sep 'm2m_' subject sep gm_mask], 0);
    elseif (tmp2~=0 && which_pipeline==3)
      gm_mask=charm_gm_mask;
      disp(['[TAP] Found ' charm_gm_mask ' and now refine alignment (the following steps may take quite a few minutes, please have patience!) ']); 
      [~, ~, ~, ~, GM_fromMesh] = load_nii_file([subjects_folder sep subject sep 'm2m_' subject sep 'surfaces' sep charm_gm_mask], 0);
    elseif (tmp2~=0 && which_pipeline==4)
      gm_mask=roast_mask_file(1).name;
      disp(['[TAP] Found ' gm_mask ' and now refine alignment (the following steps may take quite a few minutes, please have patience!) ']); 
      [~, ~, ~, ~, GM_fromMesh] = load_nii_file([subjects_folder sep subject sep gm_mask], eps);
      GM_fromMesh(find(GM_fromMesh~=2))=0;
      GM_fromMesh(find(GM_fromMesh==2))=1;
    end
    
    CC = bwconncomp(GM_fromMesh,26);
    maxi=0;
    for i=1:CC.NumObjects
        if (length(cell2mat(CC.PixelIdxList(i)))>maxi)
            brain_voxels=cell2mat(CC.PixelIdxList(i));
            maxi=length(cell2mat(CC.PixelIdxList(i)));
            ind=i;
        end
    end
    
    [gm_x gm_y gm_z] = ind2sub(size(GM_fromMesh),brain_voxels); 
    
    GM_fromMesh(:)=0;
    GM_fromMesh(brain_voxels)=1;
    gm_surf_voxel=zeros(length(gm_x),1);
    neigh=[+0, +0, +0;
    +0, +1, +0;
    +0, -1, +0;
    +1, +1, +0;
    +1, -1, +0;
    +1, +0, +0;
    -1, +1, +0;
    -1, -1, +0;
    -1, +0, +0;
    +0, +0, +1;
    +0, +1, +1;
    +0, -1, +1;
    +1, +1, +1;
    +1, -1, +1;
    +1, +0, +1;
    -1, +1, +1;
    -1, -1, +1;
    -1, +0, +1;
    +0, +0, -1;
    +0, +1, -1;
    +0, -1, -1;
    +1, +1, -1;
    +1, -1, -1;
    +1, +0, -1;
    -1, +1, -1;
    -1, -1, -1;
    -1, +0, -1];
    
    for i=1:length(gm_x)
      val=all([
          GM_fromMesh(gm_x(i)+neigh( 1,1), gm_y(i)+neigh( 1,2), gm_z(i)+neigh( 1,3)) 
          GM_fromMesh(gm_x(i)+neigh( 2,1), gm_y(i)+neigh( 2,2), gm_z(i)+neigh( 2,3))
          GM_fromMesh(gm_x(i)+neigh( 3,1), gm_y(i)+neigh( 3,2), gm_z(i)+neigh( 3,3))
          GM_fromMesh(gm_x(i)+neigh( 4,1), gm_y(i)+neigh( 4,2), gm_z(i)+neigh( 4,3))
          GM_fromMesh(gm_x(i)+neigh( 5,1), gm_y(i)+neigh( 5,2), gm_z(i)+neigh( 5,3)) 
          GM_fromMesh(gm_x(i)+neigh( 6,1), gm_y(i)+neigh( 6,2), gm_z(i)+neigh( 6,3)) 
          GM_fromMesh(gm_x(i)+neigh( 7,1), gm_y(i)+neigh( 7,2), gm_z(i)+neigh( 7,3))
          GM_fromMesh(gm_x(i)+neigh( 8,1), gm_y(i)+neigh( 8,2), gm_z(i)+neigh( 8,3))        
          GM_fromMesh(gm_x(i)+neigh( 9,1), gm_y(i)+neigh( 9,2), gm_z(i)+neigh( 9,3))         
          GM_fromMesh(gm_x(i)+neigh(10,1), gm_y(i)+neigh(10,2), gm_z(i)+neigh(10,3)) 
          GM_fromMesh(gm_x(i)+neigh(11,1), gm_y(i)+neigh(11,2), gm_z(i)+neigh(11,3))
          GM_fromMesh(gm_x(i)+neigh(12,1), gm_y(i)+neigh(12,2), gm_z(i)+neigh(12,3))
          GM_fromMesh(gm_x(i)+neigh(13,1), gm_y(i)+neigh(13,2), gm_z(i)+neigh(13,3))
          GM_fromMesh(gm_x(i)+neigh(14,1), gm_y(i)+neigh(14,2), gm_z(i)+neigh(14,3)) 
          GM_fromMesh(gm_x(i)+neigh(15,1), gm_y(i)+neigh(15,2), gm_z(i)+neigh(15,3)) 
          GM_fromMesh(gm_x(i)+neigh(16,1), gm_y(i)+neigh(16,2), gm_z(i)+neigh(16,3))
          GM_fromMesh(gm_x(i)+neigh(17,1), gm_y(i)+neigh(17,2), gm_z(i)+neigh(17,3))        
          GM_fromMesh(gm_x(i)+neigh(18,1), gm_y(i)+neigh(18,2), gm_z(i)+neigh(18,3))        
          GM_fromMesh(gm_x(i)+neigh(19,1), gm_y(i)+neigh(19,2), gm_z(i)+neigh(19,3)) 
          GM_fromMesh(gm_x(i)+neigh(20,1), gm_y(i)+neigh(20,2), gm_z(i)+neigh(20,3))
          GM_fromMesh(gm_x(i)+neigh(21,1), gm_y(i)+neigh(21,2), gm_z(i)+neigh(21,3))
          GM_fromMesh(gm_x(i)+neigh(22,1), gm_y(i)+neigh(22,2), gm_z(i)+neigh(22,3))
          GM_fromMesh(gm_x(i)+neigh(23,1), gm_y(i)+neigh(23,2), gm_z(i)+neigh(23,3)) 
          GM_fromMesh(gm_x(i)+neigh(24,1), gm_y(i)+neigh(24,2), gm_z(i)+neigh(24,3)) 
          GM_fromMesh(gm_x(i)+neigh(25,1), gm_y(i)+neigh(25,2), gm_z(i)+neigh(25,3))
          GM_fromMesh(gm_x(i)+neigh(26,1), gm_y(i)+neigh(26,2), gm_z(i)+neigh(26,3))        
          GM_fromMesh(gm_x(i)+neigh(27,1), gm_y(i)+neigh(27,2), gm_z(i)+neigh(27,3)) 
          ]);
          gm_surf_voxel(i)=not(val);
    end

    gm_surf_voxel=find(gm_surf_voxel==1);
    nr_gm_surf_voxel=length(gm_surf_voxel);
    GM_fromMesh=[gm_x(gm_surf_voxel) gm_y(gm_surf_voxel) gm_z(gm_surf_voxel)];
    GM_fromMesh=[GM_fromMesh ones(size(GM_fromMesh,1),1)]; %MRI coordinates, flip x,y
    GM_fromMesh=(GM_fromMesh*TMaskMRI);
    GM_fromMesh=GM_fromMesh(:,1:3);

    brain_nodes=brain.node';
    mod_max_y=find(brain_nodes(:,2)==max(brain_nodes(:,2)));
    mod_max_y=mean(brain_nodes(mod_max_y,:),1);
    
    T=eye(4);
    [mindis index dis] = min_distance(GM_fromMesh(:,1:3) , brain_nodes(1:fine_alignment_voxel_stepsize:end,:));
    overall_min_dis=mean(mindis);
    
    org=GM_fromMesh;
    directions=[1 -1];
    for x=1:length(directions)
        for y=1:length(directions)
            for z=1:length(directions)
               tmp=org;
               mean_tmp=mean(tmp);
               tmp(:,1)=tmp(:,1)-mean_tmp(1);
               tmp(:,2)=tmp(:,2)-mean_tmp(2);
               tmp(:,3)=tmp(:,3)-mean_tmp(3);
               tmp(:,1)=tmp(:,1)*directions(x);
               tmp(:,2)=tmp(:,2)*directions(y);
               tmp(:,3)=tmp(:,3)*directions(z);
               mean_tmp_flipped=mean(tmp);
               tmp(:,1)=tmp(:,1)+mean_tmp(1);
               tmp(:,2)=tmp(:,2)+mean_tmp(2);
               tmp(:,3)=tmp(:,3)+mean_tmp(3);
               %align y, could be also x or z+
               vox_max_y=find(tmp(:,2)==max(tmp(:,2)));
               vox_max_y=mean(tmp(vox_max_y,:),1);
               diff_y = mod_max_y-vox_max_y;
               aligned_y=tmp;
               aligned_y(:,1)=aligned_y(:,1)+diff_y(1);
               aligned_y(:,2)=aligned_y(:,2)+diff_y(2);
               aligned_y(:,3)=aligned_y(:,3)+diff_y(3);
               [mindis index dis] = min_distance(aligned_y , brain_nodes(1:alignment_voxel_stepsize:end,:));
               mindis=mean(mindis);
               if(mindis<overall_min_dis)
                  GM_fromMesh=aligned_y;
                  overall_min_dis=mindis;
               end
            end
        end
    end
    
    offset_x=0;
    offset_y=0;
    offset_z=0;
    
    for i=[-2 -1 0 1 2]
       for l=[-2 -1 0 1 2]
         for k=[-2 -1 0 1 2] 
             tmp_voxel=GM_fromMesh;
             tmp_voxel(:,1)=tmp_voxel(:,1)+i;
             tmp_voxel(:,2)=tmp_voxel(:,2)+l;
             tmp_voxel(:,3)=tmp_voxel(:,3)+k;
             tmp_brain=brain.node';
             [mindis index dis] = min_distance( tmp_voxel, tmp_brain(1:alignment_voxel_stepsize:end,:));
             mean_dis=mean(mindis);
             if (mean_dis<overall_min_dis)
                overall_min_dis=mean_dis;
                offset_x=i;
                offset_y=l;
                offset_z=k;
             end
         end
       end
    end
    
%[~,I]=sortrows(tmp_brain(:,1),'ascend');

%     voxelspace_to_mesh_moving_range=-floor(overall_min_dis):0.5:floor(overall_min_dis);
%     if (exist('save_tmp','var'))
%      [MoveMeshRot,MoveMeshT] = rigid_transform_3D(GM_fromMesh, save_tmp);
%      T(1:3,1:3)=MoveMeshRot;
%      T(4,1:3)=MoveMeshT;
%      GM_fromMesh(:,4)=1;
%      GM_fromMesh=(GM_fromMesh*T);
%      GM_fromMesh=GM_fromMesh(:,1:3);               
%     end
%     
%     nvec=trinormal(brain.face',brain.node');
%     mini=intmax;
%     tmp_brain=GM_fromMesh;
%     
%     [~,I]=sortrows(tmp_brain(:,1),'ascend');
%     minx=I(1:fine_alignment_number_voxel_to_match);
%     [~,I]=sortrows(tmp_brain(:,1),'descend');
%     maxx=I(1:fine_alignment_number_voxel_to_match);
%     [~,I]=sortrows(tmp_brain(:,2),'ascend');
%     miny=I(1:fine_alignment_number_voxel_to_match);
%     [~,I]=sortrows(tmp_brain(:,2),'descend');
%     maxy=I(1:fine_alignment_number_voxel_to_match);
%     [~,I]=sortrows(tmp_brain(:,3),'descend');
%     maxz=I(1:fine_alignment_number_voxel_to_match);
%     
%     ind=[minx; maxx; miny; maxy; maxz];
%     tmp_brain=tmp_brain(ind,:);
    
%     r = get_mesh_elm_centers(brain);
%     [mindis index dis] = min_distance(r, tmp_brain);
%     relevant_triangles=[];
%     for i=1:length(index)
%         tmp=r;
%         tmp(:,1)=tmp(:,1)-r(index(i),1);
%         tmp(:,2)=tmp(:,2)-r(index(i),2);
%         tmp(:,3)=tmp(:,3)-r(index(i),3);
%         tmp=sqrt(sum(tmp.^2,2));
%         relevant_triangles=unique([relevant_triangles; find(tmp<tri_neighborhood)]);
%     end
%     tic;
%     [Dorg,P,F] = distance_to_surfaces_vectorized(brain.face',brain.node',tmp_brain,nvec);    
%     d_org_x=[Dorg(1:2*fine_alignment_number_voxel_to_match)];
%     d_org_y=[Dorg(2*fine_alignment_number_voxel_to_match+1:4*fine_alignment_number_voxel_to_match)];
%     d_org_z=[Dorg(4*fine_alignment_number_voxel_to_match+1:5*fine_alignment_number_voxel_to_match)]; 
%     mean_org_x=abs(mean(d_org_x));
%     mean_org_y=abs(mean(d_org_y));
%     mean_org_z=abs(mean(d_org_z));  
%     mini=mean(Dorg);
%     offset_x=0; offset_y=0; offset_z=0;
%     count=1;
%     for i=voxelspace_to_mesh_moving_range
%        for l=voxelspace_to_mesh_moving_range
%          for k=voxelspace_to_mesh_moving_range
%             if (~(i==0 && l==0 && k==0))
%              tmp_voxel_moved=tmp_brain;
%              tmp_voxel_moved(:,1)=tmp_voxel_moved(:,1)+i;
%              tmp_voxel_moved(:,2)=tmp_voxel_moved(:,2)+l;
%              tmp_voxel_moved(:,3)=tmp_voxel_moved(:,3)+k;
%              [D,P,F] = distance_to_surfaces_vectorized(brain.face',brain.node',tmp_voxel_moved,nvec);           
%              d_x=[D(1:2*fine_alignment_number_voxel_to_match)];
%              d_y=[D(2*fine_alignment_number_voxel_to_match+1:4*fine_alignment_number_voxel_to_match)];
%              d_z=[D(4*fine_alignment_number_voxel_to_match+1:5*fine_alignment_number_voxel_to_match)];    
%              mean_x=abs(mean(d_x));
%              mean_y=abs(mean(d_y));
%              mean_z=abs(mean(d_z));
%              mean_dis=max([mean_x mean_y mean_z]);
%              mean_dis=mean(D);
%             % if (mean_x<=mean_org_x && mean_y<=mean_org_y && mean_z<=mean_org_z)
%                if (abs(mean_dis)<abs(mini))
%                 mini=mean_dis;
%                 offset_x=i;
%                 offset_y=l;
%                 offset_z=k;
%                end
%              %end
%             end
%          end
%        end
%     end
   
   GM_fromMesh(:,1)=GM_fromMesh(:,1)+offset_x;
   GM_fromMesh(:,2)=GM_fromMesh(:,2)+offset_y;
   GM_fromMesh(:,3)=GM_fromMesh(:,3)+offset_z; 
   
   [MoveMeshRot,MoveMeshT] = rigid_transform_3D(org, GM_fromMesh);
   T(1:3,1:3)=MoveMeshRot;
   T(4,1:3)=MoveMeshT;
   
else
    disp(['[TAP] Count not find file ' gm_mask ' or ' charm_gm_mask '. This is required. Abort.']); 
    return ;
end

%Field1.node=GM_fromMesh';
%Field1.field=zeros(1,length(GM_fromMesh));
%save([subject '_alignment.mat'],'-V6','brain','Field1');


ROI_size = abs(TMaskMRI(1:3,1:3));
ROI_size = sqrt(sum((ROI_size(find(ROI_size~=0))).^2))/2; %radius for one voxel
if ( (size(p,1)>1 && size(p,2)==3) || (size(p,2)>1 && size(p,1)==3) )
  avr_roi_voxel=mean(p);
  q=p;
  q(:,1)=q(:,1)-avr_roi_voxel(1);
  q(:,2)=q(:,2)-avr_roi_voxel(2);
  q(:,3)=q(:,3)-avr_roi_voxel(3);
  q=sqrt(sum(q.^2,2));
  ROI_size = mean(q);
end

p(:,4)=1;
simnibsModel=(p*TMaskMRI*T);
simnibsModel=simnibsModel(:,1:3);
simnibsModel(:,1)=simnibsModel(:,1)+offset_x;
simnibsModel(:,2)=simnibsModel(:,2)+offset_y;
simnibsModel(:,3)=simnibsModel(:,3)+offset_z;

if ( (size(simnibsModel,1)==1 && size(simnibsModel,2)==3) || (size(simnibsModel,1)==3 && size(simnibsModel,2)==1))
 coord = simnibsModel;
else
 coord = mean(simnibsModel);
end

if (which_pipeline==1 || which_pipeline==2 || which_pipeline==4)
  disp(['[TAP] target: ' char(mask) ' has NIfTI:aligned coordinates = [' num2str(coord(1)) ' , ' num2str(coord(2)) ' , ' num2str(coord(3)) '] in Brainsight.']);
end


