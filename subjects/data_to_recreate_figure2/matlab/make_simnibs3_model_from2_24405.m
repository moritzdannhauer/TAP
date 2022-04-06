clear all;
ANISOTROPIC=1;
addpath('/Users/rothermichk17/Applications/SimNIBS/matlab');
load /Users/rothermichk17/Pipeline/24405/24405_001_prepared_tensor.mat;
tensors=tet_mesh_ef.field';
msh=mesh_load_gmsh4('/Users/rothermichk17/S262_bin.msh');

msh.nodes=double(tet_mesh.node');
msh.tetrahedra=int32(tet_mesh.cell');
msh.tetrahedron_regions=int32(tet_mesh.field');
msh.triangles=int32(tri_mesh.face');
msh.triangle_regions=int32(tri_mesh.field');
msh.element_data{1}.tridata=double(tri_mesh_ef.field'); %this will not be used in simnibs3

%add isotropic conductivities
f=tet_mesh.field;
conductivity=double(zeros(length(f),3)); %add simnibs standard conductivities to the anisotropic brain tensors
conductivity(find(f==1),[1 5 9])=0.126; %white brain matter
conductivity(find(f==2),[1 5 9])=0.275; %gray matter
conductivity(find(f==3),[1 5 9])=1.654; %CSF
conductivity(find(f==4),[1 5 9])=0.01; %skull
conductivity(find(f==5),[1 5 9])=0.465; %skin
conductivity(find(f==6),[1 5 9])=0.5; %eye balls: mostly not used
if (ANISOTROPIC==1)
 brain_ind=find(f<=2);
 brain_tensor_avail=find(tensors(brain_ind,2)~=0);
 conductivity(brain_ind(brain_tensor_avail),:)=tensors(brain_ind(brain_tensor_avail),:);
end
msh.element_data{1}.tetdata=double(conductivity);
mesh_save_gmsh4(msh, '../24405/24405_bin.msh');
tet_mesh_ef.field=conductivity';
save -V6 ../24405/24405_bin.mat tet_mesh_ef;







