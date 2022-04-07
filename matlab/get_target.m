function [coord, which_pipeline] = get_target(subjects_folder, subject, mask, mri, sep, brain)
epsi=1e-8;
which_pipeline=0;
alignment_voxel_stepsize=1000; %use every 1000th voxel of gm mask to match brain surface

if ~exist(([subjects_folder sep subject sep mask]),'file')
  error(['[TAP] ERROR: File '  [subject sep mask] ' does not exist. ' ]);
end

if exist(([subjects_folder sep subject sep 'm2m_' subject sep 'headreco_log.html']),'file')
  disp(['[TAP] ' subject ' processed with headreco pipeline. ' ]);
  which_pipeline=1;
elseif exist(([subjects_folder sep subject sep 'm2m_' subject sep 'mri2mesh_log.html']),'file')
   if exist(([subjects_folder sep subject sep 'm2m_' subject sep 'gm_only.nii.gz']),'file')
       which_pipeline=3;
       disp(['[TAP] ' subject ' processed with mri2mesh/SimNIBS 2.x pipeline. ' ]);
   elseif exist(([subjects_folder sep subject sep 'm2m_' subject sep 'gm_fromMesh.nii.gz']),'file')
       which_pipeline=2;
       disp(['[TAP] ' subject ' processed with mri2mesh/SimNIBS 3.x pipeline. ' ]);
   else
     error(['[TAP] ERROR: Could not detect with what SimNIBS version' [subject] ' was processed. ' ]);  
   end
else
  error(['[TAP] ERROR: Could not detect with what SimNIBS version' [subject] ' was processed. ' ]);     
end

%nii files can be saved in very differently, load_untouch_nii seems to be
%most robust
Mask_MRI_nii=load_untouch_nii([subjects_folder sep subject sep mask]);

if (length(size(Mask_MRI_nii.img))~=3)
   error('[TAP] ERROR: The NII file does not contain a 3-D volume/image.'); 
end

[x y z] = ind2sub(size(Mask_MRI_nii.img),find(Mask_MRI_nii.img~=0)); 
p=[x y z]; %MRI coordinates, flip x,y

if (isempty(x) || isempty(y) || isempty(z))
  error(['[TAP] ERROR: The volume in NII file '  [subject sep mask] ' is empty, so no non-zero value could be found and used as TMS target. ' ]);
end

if (any(p<0))
    error(['[TAP] ERROR: Voxel coordinates are negative, not possible.' ]); 
end

TMaskMRI=zeros(4,4);
TMaskMRI(1,:)=Mask_MRI_nii.hdr.hist.srow_x;
TMaskMRI(2,:)=Mask_MRI_nii.hdr.hist.srow_y;
TMaskMRI(3,:)=Mask_MRI_nii.hdr.hist.srow_z;
TMaskMRI(4,4)=1;

if (length(Mask_MRI_nii.hdr.dime.dim)>4) 
 dimMaskMRI = size(Mask_MRI_nii.img);
end

if (isempty(dimMaskMRI) || any(dimMaskMRI==0))
   error(['[TAP] ERROR: NII file '  [subject sep mask] ' has no 3-D volume. ' ]); 
end

if  (any(all(TMaskMRI(1:3,1:3)==0)==1) || (Mask_MRI_nii.hdr.hist.qform_code~=0)) %this can happen, lets try something else
 if (length(Mask_MRI_nii.hdr.dime.pixdim)>4)
  if (~isempty(Mask_MRI_nii.hdr.hist.qoffset_x) && ~isempty(Mask_MRI_nii.hdr.hist.qoffset_y) && ~isempty(Mask_MRI_nii.hdr.hist.qoffset_z))
   %qform or sform
   if (Mask_MRI_nii.hdr.hist.qform_code>0) %code from niftiimage.m: line 477-504
     b = Mask_MRI_nii.hdr.hist.quatern_b;
     c = Mask_MRI_nii.hdr.hist.quatern_c;
     d = Mask_MRI_nii.hdr.hist.quatern_d;
     if 1.0-(b*b+c*c+d*d) < 0
        a = 0;
     else
        a = sqrt(1.0-(b*b+c*c+d*d));
     end
     qfactor = Mask_MRI_nii.hdr.dime.pixdim(1);
     if qfactor == 0
        qfactor = 1; 
     end
     i = Mask_MRI_nii.hdr.dime.pixdim(2);
     j = Mask_MRI_nii.hdr.dime.pixdim(3);
     k = qfactor * Mask_MRI_nii.hdr.dime.pixdim(4);
     R = [a*a+b*b-c*c-d*d 2*b*c-2*a*d     2*b*d+2*a*c
          2*b*c+2*a*d     a*a+c*c-b*b-d*d 2*c*d-2*a*b
          2*b*d-2*a*c     2*c*d+2*a*b     a*a+d*d-c*c-b*b];
     T = [Mask_MRI_nii.hdr.hist.qoffset_x, ...
                         Mask_MRI_nii.hdr.hist.qoffset_y, ...
                         Mask_MRI_nii.hdr.hist.qoffset_z];
     R = R * diag([i j k]);
     TMaskMRI=[R zeros(3,1); T 1]';
   else 
    TMaskMRI=zeros(4,4); 
    TMaskMRI(1,1) = Mask_MRI_nii.hdr.dime.pixdim(2);
    TMaskMRI(2,2) = Mask_MRI_nii.hdr.dime.pixdim(3);
    TMaskMRI(3,3) = Mask_MRI_nii.hdr.dime.pixdim(4);
    TMaskMRI(4,1) = Mask_MRI_nii.hdr.hist.qoffset_x;
    TMaskMRI(4,2) = Mask_MRI_nii.hdr.hist.qoffset_y;
    TMaskMRI(4,3) = Mask_MRI_nii.hdr.hist.qoffset_z; 
    TMaskMRI(4,4)=1;
   end
  end
 else
   error(['[TAP] ERROR: File '  [subject sep mask] ' has no image dimensions/qoffset, probably corrupt. ' ]);  
 end
end

%check if the coordinate offset is at the right place
if (any(TMaskMRI(4,1:3)==0) && any(TMaskMRI(1:3,4)~=0))
 TMaskMRI(4,1)=TMaskMRI(1,4);
 TMaskMRI(4,2)=TMaskMRI(2,4);
 TMaskMRI(4,3)=TMaskMRI(3,4);
 TMaskMRI(1,4)=0;
 TMaskMRI(2,4)=0;
 TMaskMRI(3,4)=0;  
end

Mask_MRI_CS=[]; %coordinate system convention
for i=1:3    
 ind=find(abs(TMaskMRI(i,1:3))==max(abs(TMaskMRI(i,1:3))));
 if (ind==1)
   if (TMaskMRI(i,ind)>0)
     Mask_MRI_CS=[Mask_MRI_CS 'L'];  
   else
     Mask_MRI_CS=[Mask_MRI_CS 'R'];  
   end
 end
 
 if (ind==2)
    if (TMaskMRI(i,ind)>0)
      Mask_MRI_CS=[Mask_MRI_CS 'P'];
    else
      Mask_MRI_CS=[Mask_MRI_CS 'A'];  
    end
 end
 
 if (ind==3)
    if (TMaskMRI(i,ind)>0)
      Mask_MRI_CS=[Mask_MRI_CS 'I'];  
    else
      Mask_MRI_CS=[Mask_MRI_CS 'S'];  
    end
 end
end

disp(['[TAP] Loading ' mask ' (' Mask_MRI_CS ')' ]);

if ~exist(([subjects_folder sep subject sep mri]),'file')
  error(['[TAP] ERROR: File '  [subject sep mri] ' does not exist. ' ]);
end

SimNIBS_MRI_nii = load_untouch_nii([subjects_folder sep subject sep mri]);

TsimNIBSMRI=zeros(4,4);
TsimNIBSMRI(1,:)=SimNIBS_MRI_nii.hdr.hist.srow_x;
TsimNIBSMRI(2,:)=SimNIBS_MRI_nii.hdr.hist.srow_y;
TsimNIBSMRI(3,:)=SimNIBS_MRI_nii.hdr.hist.srow_z;
TsimNIBSMRI(4,4)=1;

if (length(SimNIBS_MRI_nii.hdr.dime.dim)>4) 
 dimSimnibsMRI = size(SimNIBS_MRI_nii.img);
end

if (isempty(dimSimnibsMRI) || any(dimSimnibsMRI==0))
   error(['[TAP] ERROR: NII file '  [subject sep mri] ' has no 3-D volume. ' ]); 
end


if  (any(all(TsimNIBSMRI(1:3,1:3)==0)==1) || (SimNIBS_MRI_nii.hdr.hist.qform_code~=0)) %this can happen, lets try something else
 if (length(SimNIBS_MRI_nii.hdr.dime.pixdim)>4)
  if (~isempty(SimNIBS_MRI_nii.hdr.hist.qoffset_x) && ~isempty(SimNIBS_MRI_nii.hdr.hist.qoffset_y) && ~isempty(SimNIBS_MRI_nii.hdr.hist.qoffset_z))
   %qform or sform
   if (SimNIBS_MRI_nii.hdr.hist.qform_code>0) %code from niftiimage.m: line 477-504
     b = SimNIBS_MRI_nii.hdr.hist.quatern_b;
     c = SimNIBS_MRI_nii.hdr.hist.quatern_c;
     d = SimNIBS_MRI_nii.hdr.hist.quatern_d;
     if 1.0-(b*b+c*c+d*d) < 0
        a = 0;
     else
        a = sqrt(1.0-(b*b+c*c+d*d));
     end
     qfactor = SimNIBS_MRI_nii.hdr.dime.pixdim(1);
     if qfactor == 0
        qfactor = 1; 
     end
     i = SimNIBS_MRI_nii.hdr.dime.pixdim(2);
     j = SimNIBS_MRI_nii.hdr.dime.pixdim(3);
     k = qfactor * SimNIBS_MRI_nii.hdr.dime.pixdim(4);
     R = [a*a+b*b-c*c-d*d 2*b*c-2*a*d     2*b*d+2*a*c
          2*b*c+2*a*d     a*a+c*c-b*b-d*d 2*c*d-2*a*b
          2*b*d-2*a*c     2*c*d+2*a*b     a*a+d*d-c*c-b*b];
     T = [SimNIBS_MRI_nii.hdr.hist.qoffset_x, ...
                         SimNIBS_MRI_nii.hdr.hist.qoffset_y, ...
                         SimNIBS_MRI_nii.hdr.hist.qoffset_z];
     R = R * diag([i j k]);
     TsimNIBSMRI=[R zeros(3,1); T 1]';
   else
    TsimNIBSMRI=zeros(4,4);
    TsimNIBSMRI(1,1) = SimNIBS_MRI_nii.hdr.dime.pixdim(2);
    TsimNIBSMRI(2,2) = SimNIBS_MRI_nii.hdr.dime.pixdim(3);
    TsimNIBSMRI(3,3) = SimNIBS_MRI_nii.hdr.dime.pixdim(4);
    TsimNIBSMRI(4,1) = SimNIBS_MRI_nii.hdr.hist.qoffset_x;
    TsimNIBSMRI(4,2) = SimNIBS_MRI_nii.hdr.hist.qoffset_y;
    TsimNIBSMRI(4,3) = SimNIBS_MRI_nii.hdr.hist.qoffset_z; 
    TsimNIBSMRI(4,4)=1;
   end 
   else
    error(['[TAP] ERROR: File '  [subject sep mask] ' has no image dimensions/qoffset, probably corrupt. ' ]);  
  end
 end
end

if (any(TsimNIBSMRI(4,1:3)==0) && any(TsimNIBSMRI(1:3,4)~=0))
 TsimNIBSMRI(4,1)=TsimNIBSMRI(1,4);
 TsimNIBSMRI(4,2)=TsimNIBSMRI(2,4);
 TsimNIBSMRI(4,3)=TsimNIBSMRI(3,4);
 TsimNIBSMRI(1,4)=0;
 TsimNIBSMRI(2,4)=0;
 TsimNIBSMRI(3,4)=0;
end

SimNIBS_MRI_CS=[]; %coordinate system convention
for i=1:3    
 ind=find(abs(TsimNIBSMRI(i,1:3))==max(abs(TsimNIBSMRI(i,1:3))));
 if (ind==1)
   if (TsimNIBSMRI(i,ind)>0)
     SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'L'];  
   else
     SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'R'];  
   end
 end
 
 if (ind==2)
    if (TsimNIBSMRI(i,ind)>0)
      SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'P'];
    else
      SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'A'];  
    end
 end
 
 if (ind==3)
    if (TsimNIBSMRI(i,ind)>0)
      SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'I'];  
    else
      SimNIBS_MRI_CS=[SimNIBS_MRI_CS 'S'];  
    end
 end
end

disp(['[TAP] Loading ' mri  ' (' SimNIBS_MRI_CS ')' ]);

% if (~strcmp(Mask_MRI_CS,SimNIBS_MRI_CS))
%   disp(['[TAP] We recommand to use ' mri ' directly to transform/label the brain TMS target to/in.']);
%   error(['[TAP] Error: Please use same Nifti coordinate system definition for ' mask ' as used in ' mri ', which is (' SimNIBS_MRI_CS ')'  ]);  
% end

%fine alignment MRI and mesh
if (which_pipeline==1 || which_pipeline==2)
  if (exist([subjects_folder sep subject sep 'm2m_' subject sep 'gm_fromMesh.nii.gz'],'file')~=0) 
    disp(['[TAP] Found gm_fromMesh.nii.gz and now refine alignment (this may take a few minutes) ']);
    GM_fromMesh=load_untouch_nii([subjects_folder sep subject sep 'm2m_' subject sep 'gm_fromMesh.nii.gz']);
    GM_fromMesh=GM_fromMesh.img;
    [gm_x gm_y gm_z] = ind2sub(size(GM_fromMesh),find(GM_fromMesh~=0)); 
    gm_surf_voxel=zeros(length(gm_x),1);
    for i=1:length(gm_x)
      val=all([
          GM_fromMesh(gm_x(i)+0, gm_y(i)+0, gm_z(i)+0) 
          GM_fromMesh(gm_x(i)+0, gm_y(i)+1, gm_z(i)+0)
          GM_fromMesh(gm_x(i)+0, gm_y(i)-1, gm_z(i)+0)
          GM_fromMesh(gm_x(i)+1, gm_y(i)+1, gm_z(i)+0)
          GM_fromMesh(gm_x(i)+1, gm_y(i)-1, gm_z(i)+0) 
          GM_fromMesh(gm_x(i)+1, gm_y(i)+0, gm_z(i)+0) 
          GM_fromMesh(gm_x(i)-1, gm_y(i)+1, gm_z(i)+0)
          GM_fromMesh(gm_x(i)-1, gm_y(i)-1, gm_z(i)+0)        
          GM_fromMesh(gm_x(i)-1, gm_y(i)+0, gm_z(i)+0)         
          GM_fromMesh(gm_x(i)+0, gm_y(i)+0, gm_z(i)+1) 
          GM_fromMesh(gm_x(i)+0, gm_y(i)+1, gm_z(i)+1)
          GM_fromMesh(gm_x(i)+0, gm_y(i)-1, gm_z(i)+1)
          GM_fromMesh(gm_x(i)+1, gm_y(i)+1, gm_z(i)+1)
          GM_fromMesh(gm_x(i)+1, gm_y(i)-1, gm_z(i)+1) 
          GM_fromMesh(gm_x(i)+1, gm_y(i)+0, gm_z(i)+1) 
          GM_fromMesh(gm_x(i)-1, gm_y(i)+1, gm_z(i)+1)
          GM_fromMesh(gm_x(i)-1, gm_y(i)-1, gm_z(i)+1)        
          GM_fromMesh(gm_x(i)-1, gm_y(i)+0, gm_z(i)+1)        
          GM_fromMesh(gm_x(i)+0, gm_y(i)+0, gm_z(i)-1) 
          GM_fromMesh(gm_x(i)+0, gm_y(i)+1, gm_z(i)-1)
          GM_fromMesh(gm_x(i)+0, gm_y(i)-1, gm_z(i)-1)
          GM_fromMesh(gm_x(i)+1, gm_y(i)+1, gm_z(i)-1)
          GM_fromMesh(gm_x(i)+1, gm_y(i)-1, gm_z(i)-1) 
          GM_fromMesh(gm_x(i)+1, gm_y(i)+0, gm_z(i)-1) 
          GM_fromMesh(gm_x(i)-1, gm_y(i)+1, gm_z(i)-1)
          GM_fromMesh(gm_x(i)-1, gm_y(i)-1, gm_z(i)-1)        
          GM_fromMesh(gm_x(i)-1, gm_y(i)+0, gm_z(i)-1) 
          ]);
          gm_surf_voxel(i)=not(val);
    end
    gm_surf_voxel=find(gm_surf_voxel==1);
    nr_gm_surf_voxel=length(gm_surf_voxel);
    GM_fromMesh=[gm_x(gm_surf_voxel) gm_y(gm_surf_voxel) gm_z(gm_surf_voxel)];
    GM_fromMesh=[GM_fromMesh ones(size(GM_fromMesh,1),1)]; %MRI coordinates, flip x,y
    GM_fromMesh=(GM_fromMesh*TMaskMRI);
    GM_fromMesh=GM_fromMesh(:,1:3);
    tic;
    mini=1e300;
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
             if (mean_dis<mini)
                mini=mean_dis;
                offset_x=i;
                offset_y=l;
                offset_z=k;
             end
         end
       end
    end
    toc
   GM_fromMesh(:,1)=GM_fromMesh(:,1)+offset_x;
   GM_fromMesh(:,2)=GM_fromMesh(:,2)+offset_y;
   GM_fromMesh(:,3)=GM_fromMesh(:,3)+offset_z;    
  else
    disp(['[TAP] Could not find gm_fromMesh.nii.gz in subfolder m2m_' subject sep ' - that means no fine adjustment of the MRI-mesh alignment possible.']);
  end  
end


disp(['[TAP] Refined alignment of grey matter segmentation to mesh difference: ' num2str(mean_dis) '[mm]']);


p(:,4)=1;
simnibsModel=(p*TMaskMRI);
simnibsModel=simnibsModel(:,1:3);
simnibsModel(:,1)=simnibsModel(:,1)+offset_x;
simnibsModel(:,2)=simnibsModel(:,2)+offset_y;
simnibsModel(:,3)=simnibsModel(:,3)+offset_z;
        
if ( (size(simnibsModel,1)==1 && size(simnibsModel,2)==3) || (size(simnibsModel,1)==3 && size(simnibsModel,2)==1))
 coord = simnibsModel;
else
 coord = mean(simnibsModel);
end

if (which_pipeline==1 || which_pipeline==2)
  disp([' [TAP] target: ' char(mask) ' has NIfTI:aligned coordinates = [' num2str(coord(1)) ' , ' num2str(coord(2)) ' , ' num2str(coord(3)) '] in Brainsight.']);
end



