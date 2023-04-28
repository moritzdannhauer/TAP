function out = get_node_normals(input)
max_tri_neighbors=50;
face_normals=input.field'; %DEF: should point outwards
node_face_connection=zeros(length(input.node),max_tri_neighbors);
node_normals=zeros(length(input.face),3);
f=input.face';
node_face_connection(:,1)=1;
for i=1:length(input.face)
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

node_normals=zeros(length(input.node),3);
%determine mesh node normals
for i=1:length(input.node)
 tmp=face_normals(node_face_connection(i,2:2+node_face_connection(i,1)-2),:);
 node_normals(i,:)=mean(tmp); 
 node_normals(i,:)=node_normals(i,:)/norm(node_normals(i,:));
end

out=input;
out.field=node_normals';


