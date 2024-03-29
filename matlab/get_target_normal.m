function [target, scalp_normal] = get_target_normal(subjfolder, sep, subj, coord, brain, scalp, head_model_mesh, which_pipeline, brain_vert_neigh)
max_tri_neighbors=50;
max_nr_node_curv_check=20;
max_ntermeshnode_distance=5;
check_tri_neigh=2;
headreco_steps_for_crawling_along_inward_normal=20;

[dis, normal_index, ~] = min_distance(head_model_mesh.node', coord);
if (any(dis)>max_ntermeshnode_distance)
  disp('Error: projected nodes are probably far outside the mesh.');
  return ; 
end

if (which_pipeline==2)
%Idea: Looking for curvature closer to '0' to be part of a suci wall by
%loading the surface curvature generated by freesurfer 
 if (exist([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'lh.sulc'],'file') && exist([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'rh.sulc'],'file') && exist([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'lh.pial'],'file') && exist([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'rh.pial'],'file'))
  [l_curv, ~] = read_curv([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'lh.sulc']);
  [l_vert, l_f] = read_surf([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'lh.pial']);
  [r_curv, ~] = read_curv([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'rh.sulc']);
  [r_vert, r_f] = read_surf([subjfolder sep subj sep 'fs_' subj sep 'surf' sep 'rh.pial']);
  [dis1, ~, ~] = min_distance(l_vert, coord);
  [dis2, ~, ~] = min_distance(r_vert, coord);
  %find the node that is within normal_search_radius and closest to curvature
  %to 0 indicating sulci wall
  if (dis1<dis2)
    freesurfer_n=l_vert;
    freesurfer_c=l_curv;
    if (min(l_f(:))==0 && max(l_f(:))==length(freesurfer_n)-1)
      freesurfer_f=l_f+1;  
    else
      freesurfer_f=l_f;  
     end
   else
     freesurfer_n=r_vert;
     freesurfer_c=r_curv;
     if (min(r_f(:))==0 && max(r_f(:))==length(freesurfer_n)-1)
       freesurfer_f=r_f+1;  
     else
      freesurfer_f=r_f;  
     end  
  end
 
  [~, ind, ~] = min_distance(freesurfer_n, coord);
  ind1=ismember(freesurfer_f(:,1),ind);
  ind2=ismember(freesurfer_f(:,2),ind);
  ind3=ismember(freesurfer_f(:,3),ind);
  f=[ind1 + ind2 + ind3];
  f=find(f>0);
  mini_curv=1e100;
  for i=1:check_tri_neigh
    node_ind=freesurfer_f(f,:);
    node_ind=node_ind(:);
    node_ind=unique(node_ind);
    [current_c,tmp_ind]=min(abs(freesurfer_c(node_ind)));
    ind=node_ind(tmp_ind);
    if (current_c<mini_curv)
       mini_curv=current_c;
       mini_ind=ind;
    end
    ind1=ismember(freesurfer_f(:,1),node_ind);
    ind2=ismember(freesurfer_f(:,2),node_ind);
    ind3=ismember(freesurfer_f(:,3),node_ind);
    ind=ind1+ind2+ind3;
    f=find(ind>0);
  end
  
  n=brain.node';
  freesurfer_sulc_node=freesurfer_n(mini_ind,:);
  [~, normal_index, ~] = min_distance(n, freesurfer_sulc_node);
  [~, index, ~] = min_distance(n, coord);
  target.node=n(index,:)';
  node_normals=brain.field';
  target.field=node_normals(normal_index,:)';
 else
     disp('[TAP] Error: Could not find subjects freesurfer required files (rh.sulc and lh.sulc). Abort.');
 end
else
  tmp=scalp.node';
  scalp_node_normals=scalp.field';
  [~, index, ~] = min_distance(tmp, coord);
  spp_normal=scalp_node_normals(index,:);
  spp = tmp(index,:);
  spp_target_vector = coord - spp;
  spp_target_vector = spp_target_vector/norm(spp_target_vector);
  tmp = coord;
  brain_nodes = brain.node';
  brain_node_normals = brain.field'; 
  [~, coord_index, ~] = min_distance(brain_nodes, coord);
  min_angle=1e400;
  angle=abs(acosd(dot(spp_target_vector, brain_node_normals(coord_index,:)/norm(brain_node_normals(coord_index,:))))-90);
  min_angle=angle;
  save_normal=brain_node_normals(coord_index,:);
  for i=1:headreco_steps_for_crawling_along_inward_normal    
    coord_neighbors=brain_vert_neigh(coord_index,2:2+brain_vert_neigh(coord_index,1)-1);
    brain_coord_neighbors=brain_nodes(full(coord_neighbors),:);
    brain_coord_neighbors(:,1)=brain_coord_neighbors(:,1)-coord(1);
    brain_coord_neighbors(:,2)=brain_coord_neighbors(:,2)-coord(2);
    brain_coord_neighbors(:,3)=brain_coord_neighbors(:,3)-coord(3);
    for j=1:size(brain_coord_neighbors,1)
        brain_coord_neighbors(j,:)=brain_coord_neighbors(j,:)/norm(brain_coord_neighbors(j,:));
    end
    dot_res = zeros(length(brain_coord_neighbors),1);
    for j=1:length(brain_coord_neighbors)
       dot_res(j)=dot(spp_target_vector,brain_coord_neighbors(j,:));
    end
    crawl_up=find(min(dot_res)==dot_res);
    coord_index=coord_neighbors(crawl_up);
    angle=abs(acosd(dot(spp_target_vector, brain_node_normals(coord_index,:)))-90);
    if (angle < min_angle)
        min_angle=angle;
        save_normal=brain_node_normals(coord_index,:);
    end 
  end
  
  
  for i=1:headreco_steps_for_crawling_along_inward_normal    
    coord_neighbors=brain_vert_neigh(coord_index,2:2+brain_vert_neigh(coord_index,1)-1);
    brain_coord_neighbors=brain_nodes(full(coord_neighbors),:);
    brain_coord_neighbors(:,1)=brain_coord_neighbors(:,1)-coord(1);
    brain_coord_neighbors(:,2)=brain_coord_neighbors(:,2)-coord(2);
    brain_coord_neighbors(:,3)=brain_coord_neighbors(:,3)-coord(3);
    for j=1:size(brain_coord_neighbors,1)
        brain_coord_neighbors(j,:)=brain_coord_neighbors(j,:)/norm(brain_coord_neighbors(j,:));
    end
    dot_res = zeros(length(brain_coord_neighbors),1);
    for j=1:length(brain_coord_neighbors)
       dot_res(j)=dot(spp_target_vector,brain_coord_neighbors(j,:));
    end
    crawl_down=find(max(dot_res)==dot_res);
    coord_index=coord_neighbors(crawl_down);
    angle=abs(acosd(dot(spp_target_vector, brain_node_normals(coord_index,:)/norm(brain_node_normals(coord_index,:))))-90);
    if (angle < min_angle)
        min_angle=angle;
        save_normal=brain_node_normals(coord_index,:);
    end 
  end
  
  [~, coord_index, ~] = min_distance(brain_nodes, coord);
  target.node=brain_nodes(coord_index,:)';
  target.field=save_normal';
end


n=scalp.node';
f=scalp.field';
%get scalp-projected point
[~, ind, ~] = min_distance(n, coord);
n=n(ind,:);
f=f(ind,:);
f=f/norm(f);
p=n;
q=n+target.field';
qproj=q-dot(q-p,f)*f;
o=qproj-p;
o=o/norm(o);
target.field=o';

scalp_normal.node=n';
scalp_normal.field=f';
disp('Preparing surfaces and target: done.');
