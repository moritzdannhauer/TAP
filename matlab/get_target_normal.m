function [brain, scalp, target] = get_target_normal(subj, coord, material)
max_tri_neighbors=50;
%TODO add sulci wall detection

msh=mesh_load_gmsh4([subj '/' subj '_bin.msh']);
ind=find(msh.tetrahedron_regions==material);
if (isempty(ind))
  disp('Error: material label not in mesh!');
  return ;
end
[tri, nhat]=surftri(msh.nodes,msh.tetrahedra(find(msh.tetrahedron_regions==material),:));
brain.node=msh.nodes';
brain.face=tri';
brain.field=nhat';
face_normals=nhat; %DEF: should point outwards
brain = clean_tri_mesh(brain.node',brain.face',brain.field');
node_face_connection=zeros(length(brain.node),max_tri_neighbors);
node_normals=zeros(length(brain.face),3);
f=brain.face';
node_face_connection(:,1)=1;
for i=1:length(brain.face)
    t=f(i,:);
    if (any(t>0) && any(t<length(node_face_connection)))
     node_face_connection(t(1),node_face_connection(t(1),1)+1)=i;
     node_face_connection(t(1),1)=node_face_connection(t(1),1)+1;
    
     node_face_connection(t(2),node_face_connection(t(2),1)+1)=i;
     node_face_connection(t(2),1)=node_face_connection(t(2),1)+1;
    
     node_face_connection(t(3),node_face_connection(t(3),1)+1)=i;
     node_face_connection(t(3),1)=node_face_connection(t(3),1)+1;
    else
      disp('ERROR: Triangle definition corrupted!!!');
      return;
    end
end
node_face_connection=node_face_connection(:,1:max(node_face_connection(:,1)));

node_normals=zeros(length(brain.node),3);
%determine mesh node normals
% maxi=0;count=0;
for i=1:length(brain.node)
 tmp=face_normals(node_face_connection(i,2:2+node_face_connection(i,1)-2),:);
%  for j=1:length(tmp)
%    for k=1:length(tmp)  
%      a=acosd(dot(tmp(j,:),tmp(k,:)));
%      if (a>maxi)
%         maxi=a;
%         ind=i;
%         count=count+1;
%      end
%    end
%  end
 node_normals(i,:)=mean(tmp); 
 node_normals(i,:)=node_normals(i,:)/norm(node_normals(i,:));
end

[l_curv, l_fnum] = read_curv([subj '/fs_' subj '/surf/lh.sulc']);
[l_vert, l_face] = read_surf([subj '/fs_' subj '/surf/lh.pial']);
[r_curv, r_fnum] = read_curv([subj '/fs_' subj '/surf/rh.sulc']);
[r_vert, r_face] = read_surf([subj '/fs_' subj '/surf/rh.pial']);

[mindis index dis] = min_distance(brain.node', coord);
target.node=brain.node(:,index);
target.field=node_normals(index,:);


[tri, nhat]=surftri(msh.nodes,msh.tetrahedra(find(msh.tetrahedron_regions==max(msh.tetrahedron_regions)),:));
scalp.node=msh.nodes';
scalp.face=tri';
scalp.field=nhat';
disp('Preparing surfaces and target: done.');
