function p = convcoord_Voxel_TO_Nii(coord, rascorner, dims, voxel_dim)

coord=coord+1;
p=[(coord(:,1)-sign(voxel_dim(1))*rascorner(1))/abs(voxel_dim(1)) (coord(:,2)-sign(voxel_dim(2))*rascorner(2))/abs(voxel_dim(2)) (coord(:,3)-sign(voxel_dim(3))*rascorner(3))/abs(voxel_dim(3))]; %world coordinates = nifti





