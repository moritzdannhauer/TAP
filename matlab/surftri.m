function [tri, nhat, inner_node]=surftri(p,t)
%SURFTRI Find surface triangles from tetrahedra mesh
%   TRI=SURFTRI(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

nte=size(t);
% Form all faces, non-duplicates are surface triangles
faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
%node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
cell_ind=[(1:nte(1))';(1:nte(1))';(1:nte(1))';(1:nte(1))'];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
[vec, bin]=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);
%node4=node4(ix(qx));
c=t(cell_ind(ix(qx),:),:);
inner_node=cell2mat(cellfun(@setdiff, num2cell(c, 2), num2cell(tri, 2), 'UniformOutput', false));
% Orientation
v1=p(tri(:,2),:)-p(tri(:,1),:);
v2=p(tri(:,3),:)-p(tri(:,1),:);
v3=p(inner_node,:)-p(tri(:,1),:);
ix=find(dot(cross(v1,v2,2),v3,2)>0);
tri(ix,[2,3])=tri(ix,[3,2]);
v1=p(tri(:,2),:)-p(tri(:,1),:);
v2=p(tri(:,3),:)-p(tri(:,1),:);
nhat=cross(v1,v2,2);
for i=1:length(nhat)
  nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
end






