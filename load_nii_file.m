function [Transformation, p, Mask_MRI_CS, dimMaskMRI, Voxels] = load_nii_file(file_path, thres)

MRI_nii=load_untouch_nii(file_path);

if (length(size(MRI_nii.img))~=3)
   error('[TAP] ERROR: The NII file does not contain a 3-D volume/image.'); 
end

Voxels=MRI_nii.img;
if (isinf(thres))
    thres=max(Voxels(:));
end
[x y z] = ind2sub(size(Voxels),find(abs(Voxels)>=thres)); 
p=[x y z]; %MRI coordinates, flip x,y

if (isempty(x) || isempty(y) || isempty(z))
  error(['[TAP] ERROR: The volume in NII file ' file_path ' is empty, so no non-zero value could be found and used as TMS target. ' ]);
end

if (any(p<0))
    error(['[TAP] ERROR: Voxel coordinates are negative, not possible.' ]); 
end

Transformation=zeros(4,4);
Transformation(1,:)=MRI_nii.hdr.hist.srow_x;
Transformation(2,:)=MRI_nii.hdr.hist.srow_y;
Transformation(3,:)=MRI_nii.hdr.hist.srow_z;
Transformation(4,4)=1;

if (length(MRI_nii.hdr.dime.dim)>4) 
 dimMaskMRI = size(MRI_nii.img);
end

if (isempty(dimMaskMRI) || any(dimMaskMRI==0))
   error(['[TAP] ERROR: NII file '  file_path ' has no 3-D volume. ' ]); 
end

if  (any(all(Transformation(1:3,1:3)==0)==1) || (MRI_nii.hdr.hist.qform_code~=0)) %this can happen, lets try something else
 if (length(MRI_nii.hdr.dime.pixdim)>4)
  if (~isempty(MRI_nii.hdr.hist.qoffset_x) && ~isempty(MRI_nii.hdr.hist.qoffset_y) && ~isempty(MRI_nii.hdr.hist.qoffset_z))
   %qform or sform
   if (MRI_nii.hdr.hist.qform_code>0) %code from niftiimage.m: line 477-504
     b = MRI_nii.hdr.hist.quatern_b;
     c = MRI_nii.hdr.hist.quatern_c;
     d = MRI_nii.hdr.hist.quatern_d;
     if 1.0-(b*b+c*c+d*d) < 0
        a = 0;
     else
        a = sqrt(1.0-(b*b+c*c+d*d));
     end
     qfactor = MRI_nii.hdr.dime.pixdim(1);
     if qfactor == 0
        qfactor = 1; 
     end
     i = MRI_nii.hdr.dime.pixdim(2);
     j = MRI_nii.hdr.dime.pixdim(3);
     k = qfactor * MRI_nii.hdr.dime.pixdim(4);
     R = [a*a+b*b-c*c-d*d 2*b*c-2*a*d     2*b*d+2*a*c
          2*b*c+2*a*d     a*a+c*c-b*b-d*d 2*c*d-2*a*b
          2*b*d-2*a*c     2*c*d+2*a*b     a*a+d*d-c*c-b*b];
     T = [MRI_nii.hdr.hist.qoffset_x, ...
                         MRI_nii.hdr.hist.qoffset_y, ...
                         MRI_nii.hdr.hist.qoffset_z];
     R = R * diag([i j k]);
     Transformation=[R zeros(3,1); T 1]';
   else 
    Transformation=zeros(4,4); 
    Transformation(1,1) = MRI_nii.hdr.dime.pixdim(2);
    Transformation(2,2) = MRI_nii.hdr.dime.pixdim(3);
    Transformation(3,3) = MRI_nii.hdr.dime.pixdim(4);
    Transformation(4,1) = MRI_nii.hdr.hist.qoffset_x;
    Transformation(4,2) = MRI_nii.hdr.hist.qoffset_y;
    Transformation(4,3) = MRI_nii.hdr.hist.qoffset_z; 
    Transformation(4,4)=1;
   end
  end
 else
   error(['[TAP] ERROR: File '  file_path ' has no image dimensions/qoffset, probably corrupt. ' ]);  
 end
end

%check if the coordinate offset is at the right place
if (any(Transformation(4,1:3)==0) && any(Transformation(1:3,4)~=0))
 Transformation(4,1)=Transformation(1,4);
 Transformation(4,2)=Transformation(2,4);
 Transformation(4,3)=Transformation(3,4);
 Transformation(1,4)=0;
 Transformation(2,4)=0;
 Transformation(3,4)=0;  
end

Mask_MRI_CS=[]; %coordinate system convention
for i=1:3    
 ind=find(abs(Transformation(i,1:3))==max(abs(Transformation(i,1:3))));
 if (ind==1)
   if (Transformation(i,ind)>0)
     Mask_MRI_CS=[Mask_MRI_CS 'R'];  
   else
     Mask_MRI_CS=[Mask_MRI_CS 'L'];  
   end
 end
 
 if (ind==2)
    if (Transformation(i,ind)>0)
      Mask_MRI_CS=[Mask_MRI_CS 'A'];
    else
      Mask_MRI_CS=[Mask_MRI_CS 'P'];  
    end
 end
 
 if (ind==3)
    if (Transformation(i,ind)>0)
      Mask_MRI_CS=[Mask_MRI_CS 'S'];  
    else
      Mask_MRI_CS=[Mask_MRI_CS 'I'];  
    end
 end
end

disp(['[TAP] Loading ' file_path ' (' Mask_MRI_CS ')' ]);