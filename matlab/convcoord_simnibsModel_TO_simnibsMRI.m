function simnibsMRI = convcoord_simnibsModel_TO_simnibsMRI(subject, mask_name, coordinates, dimMaskMRI, dimSimnibsMRI, TMaskMRI, TsimNIBSMRI)


offset = [dimSimnibsMRI(1)/2 dimSimnibsMRI(2)/2 dimSimnibsMRI(3)/2 1];
SimNIBS_center=offset*TsimNIBSMRI;
SimNIBS_center=SimNIBS_center(1:3);

% for i=1:size(real_world,1)
%  if (SimNIBS_mesh_center(1)<real_world(i,1))
%    real_world(i,1)=SimNIBS_mesh_center(1)-real_world(i,1)+TMaskMRI(1,1); 
%  else
%    real_world(i,1)=real_world(i,1)-SimNIBS_mesh_center(1)-TMaskMRI(1,1);
%  end
%  if (SimNIBS_mesh_center(2)<real_world(i,2))
%     real_world(i,2)=SimNIBS_mesh_center(2)-real_world(i,2)+TMaskMRI(2,2);
%  else
%     real_world(i,2)=real_world(i,2)-SimNIBS_mesh_center(2)-TMaskMRI(2,2);
%  end
%  if (SimNIBS_mesh_center(3)>real_world(i,3))
%    real_world(i,3)=SimNIBS_mesh_center(3)-real_world(i,3); 
%  else
%    real_world(i,3)=real_world(i,3)-SimNIBS_mesh_center(3);  
%  end
% end

vec_from_origin=coordinates+SimNIBS_center;
vec_from_origin(:,1)=vec_from_origin(:,1)+TMaskMRI(1,1);
vec_from_origin(:,2)=vec_from_origin(:,2)+TMaskMRI(2,2);
real_world=vec_from_origin;
%real_world(:,1)=SimNIBS_mesh_center(1)-real_world(:,1)+TMaskMRI(1,1); 
%real_world(:,2)=SimNIBS_mesh_center(2)-real_world(:,2)+TMaskMRI(2,2); 
%real_world(:,3)=real_world(:,3)+SimNIBS_mesh_center(3); 

real_world(:,4)=1;
try
  real_world=real_world*inv(TMaskMRI);
 catch
  disp(['[TAP] ERROR: Transfermatrix converting targeting mask coordinates back to MRI coordinates not invertable.']) 
end


real_world=real_world(:,1:3);

real_world=abs(real_world); %mri coordinates are positive
simnibsMRI=real_world;










