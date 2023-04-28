function p = convcoord_Voxel_TO_Nii(coord, rascorner, dims, voxel_dim)

p=[coord(:,1) coord(:,2) coord(:,3)]; 
p=[p(:,1)*abs(voxel_dim(1))-rascorner(1) p(:,2)*abs(voxel_dim(2))+rascorner(2) p(:,3)*abs(voxel_dim(3))+rascorner(3)]-1; %world coordinates = nifti




