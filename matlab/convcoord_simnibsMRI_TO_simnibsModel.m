function simnibsModel = simnibsMRI2simnibsModel(subject, mask, coordinates, dimMaskMRI, dimSimnibsMRI, TMaskMRI, TsimNIBSMRI)

MRI_center = dimMaskMRI ./ 2;
SimNIBS_center = dimSimnibsMRI ./ 2;
MRI_center_neigh = [];
for i=MRI_center(1)-1:MRI_center(1)+1
  for j=MRI_center(2)-1:MRI_center(2)+1
    for k=MRI_center(3)-1:MRI_center(3)+1
        MRI_center_neigh=[MRI_center_neigh; i j k];
    end
  end
end
MRI_center_neigh = [MRI_center; MRI_center_neigh];
SimNIBS_center_neigh = [];
for i=SimNIBS_center(1)-1:SimNIBS_center(1)+1
  for j=SimNIBS_center(2)-1:SimNIBS_center(2)+1
    for k=SimNIBS_center(3)-1:SimNIBS_center(3)+1
        SimNIBS_center_neigh=[SimNIBS_center_neigh; i j k];
    end
  end
end
SimNIBS_center_neigh = [SimNIBS_center; SimNIBS_center_neigh];

SimNIBS_center_neigh(:,4)=1;
MRI_center_neigh(:,4)=1;

MRI_center_axis=(TMaskMRI*MRI_center_neigh')';
MRI_center_axis=MRI_center_axis(:,1:3);
SimNIBS_center_axis=(TsimNIBSMRI*SimNIBS_center_neigh')';
SimNIBS_center_axis=SimNIBS_center_axis(:,1:3);

[R,t] = rigid_transform_3D(MRI_center_axis, SimNIBS_center_axis);
T_Mask_to_SimNIBS = zeros(4,4);
T_Mask_to_SimNIBS(1:3,1:3) = R;
T_Mask_to_SimNIBS(1:3,4) = t';
T_Mask_to_SimNIBS(4,4) = 1;

real_world_coordinates=coordinates;
real_world_coordinates(:,4)=1;
real_world_coordinates=(TMaskMRI*real_world_coordinates')';
real_world_coordinates=real_world_coordinates(:,1:3);

%TODO: this wont work if transformation of axis span across multiple axis, i.e.,
%matrix entries are not '1 0 0' rather '0.7071 0.7071 0' 
displacement_correction=zeros(3,1);  
for i=1:3
  displacement_correction(i,1)=TMaskMRI(i,find(abs(TMaskMRI(i,1:3))==max(abs(TMaskMRI(i,1:3)))));
end
for i=1:3
   real_world_coordinates(:,i)=real_world_coordinates(:,i)-displacement_correction(i); 
end

%and now consider SimNIBS coordinate transformation
real_world_coordinates(:,4)=1;
real_world_coordinates=(T_Mask_to_SimNIBS*real_world_coordinates')';
real_world_coordinates=real_world_coordinates(:,1:3);

simnibsModel=real_world_coordinates;

%%%OLD%%%%

% % real_world=coordinates;
% % real_world(:,4)=1;
% % real_world=real_world*TMaskMRI;
% % real_world=real_world(:,1:3);
% % 
% % offset = [dimSimnibsMRI(1)/2 dimSimnibsMRI(2)/2 dimSimnibsMRI(3)/2 1];
% % SimNIBS_mesh_center=offset*TsimNIBSMRI;
% % SimNIBS_mesh_center=SimNIBS_mesh_center(1:3);

% figure; hold on;
% str={'RPA','LPA','IN','NA','CZ'};
% for i=1:length(real_world)
%  x=real_world(i,1);
%  y=real_world(i,2);
%  z=real_world(i,3);
%  plot3(x,y,z,'b.');
%  text(x,y,z,char(str{i}));
% end
% x=SimNIBS_mesh_center(1);
% y=SimNIBS_mesh_center(2);
% z=SimNIBS_mesh_center(3);
% plot3(x,y,z,'ro');
% text(x,y,z,'Center');
% axis equal;
% grid on;
% xlabel('x');
% ylabel('y');
% zlabel('z');

%for i=1:length(real_world)
 %determine case
% %  vec_from_origin=real_world-SimNIBS_mesh_center;
% %  vec_from_origin(:,1)=vec_from_origin(:,1)-TMaskMRI(1,1);
% %  vec_from_origin(:,2)=vec_from_origin(:,2)-TMaskMRI(2,2);
%  what_quadrant=sign(vec_from_origin);
%  if all(what_quadrant==[1 -1 -1]) %case 8
%       real_world(i,1)=real_world(i,1)-SimNIBS_mesh_center(1);
%       real_world(i,2)=real_world(i,2)-SimNIBS_mesh_center(2);
%   end
 
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
% real_world(i,:) = vec_from_origin(i,:);
%end

% % simnibsModel=vec_from_origin;








