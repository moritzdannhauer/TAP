function [brain, scalp, head_model_mesh, brain_vert_neigh] = get_meshSurfs(subjects_folder, subj, material, sep)

smooth_normals_brain=8;
smooth_normals_brain_weight=0.7;
smooth_normals_brain_maxdisplacement=1.0;
smooth_normals_scalp=8;
smooth_normals_scalp_weight=0.7;
smooth_normals_scalp_maxdisplacement=1.0;

msh=mesh_load_gmsh4([subjects_folder sep subj sep subj '.msh']);

ind=find(msh.tetrahedron_regions==material);
if (isempty(ind))
  disp('[TAP] Error: material label not in mesh!');
  return ;
end

head_model_mesh.node=msh.nodes';
if(size(head_model_mesh.node,1)~=3 &&  size(head_model_mesh.node,2)==3)
    head_model_mesh.node=head_model_mesh.node';
end
head_model_mesh.cell=msh.tetrahedra';
if(size(head_model_mesh.cell,1)~=4 &&  size(head_model_mesh.cell,2)==4)
    head_model_mesh.cell=head_model_mesh.cell';
end
head_model_mesh.field=msh.tetrahedron_regions';
if(size(head_model_mesh.field,1)~=1 &&  size(head_model_mesh.field,2)==1)
    head_model_mesh.field=head_model_mesh.field';
end
if (size(head_model_mesh.cell,2)~=size(head_model_mesh.field,2))
  disp('[TAP] Error: mesh elements and labels are not equal. aborting TAP !');  
end

[tri, nhat, inner_node]=surftri(head_model_mesh.node',head_model_mesh.cell(:,ismember(head_model_mesh.field',material))');
all_nodes=msh.nodes;
brain.node=all_nodes';
brain.face=tri';
brain.field=inner_node';
brain = clean_tri_mesh(brain.node',brain.face',brain.field');

%Smooth Brain to get a better estimate of the normal = sulci wall
[vert_neigh vert_tri_neigh] = neighborelem(brain.node',brain.face');
brain_vert_neigh = vert_neigh;
n=brain.node';
f=brain.face';
%check if normals are really pointing outwards
nvec=-trinormal(brain.face',brain.node');
for i=1:length(brain.face)
    n1=mean(n(f(i,:),:))-all_nodes(inner_node(i),:);
    n1=n1/norm(n1);
    if (sign(dot(n1,nvec(i,:)))<0)
        nvec(i,:)=-nvec(i,:);
    end
end
brain.field=nvec';
brain=get_node_normals(brain);
[newVertex, newNormals] = SmoothLaplace(sparse(brain.node'), vert_neigh, brain.field', smooth_normals_brain, smooth_normals_brain_weight, [], smooth_normals_brain_maxdisplacement);
for i=1:size(newNormals,1)
   newNormals(i,:)=newNormals(i,:)/norm(newNormals(i,:)); 
end
brain.field=newNormals';

[tri, nhat, inner_node]=surftri(head_model_mesh.node',head_model_mesh.cell');
all_nodes=msh.nodes;
scalp.node=all_nodes';
scalp.face=tri';
scalp.field=inner_node';
scalp = clean_tri_mesh(scalp.node',scalp.face',scalp.field');

%Smooth Brain to get a better estimate of the normal = sulci wall
[vert_neigh vert_tri_neigh] = neighborelem(scalp.node',scalp.face');
n=scalp.node';
f=scalp.face';
%check if normals are really pointing outwards
nvec=-trinormal(scalp.face',scalp.node');
for i=1:length(scalp.face)
    n1=mean(n(f(i,:),:))-all_nodes(inner_node(i),:);
    n1=n1/norm(n1);
    if (sign(dot(n1,nvec(i,:)))<0)
        nvec(i,:)=-nvec(i,:);
    end
end
scalp.field=nvec';
scalp=get_node_normals(scalp);
[newVertex, newNormals] = SmoothLaplace(sparse(scalp.node'), vert_neigh, scalp.field', smooth_normals_brain, smooth_normals_brain_weight, [], smooth_normals_brain_maxdisplacement);
for i=1:size(newNormals,1)
   newNormals(i,:)=newNormals(i,:)/norm(newNormals(i,:)); 
end
scalp.field=newNormals';






