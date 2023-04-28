function [newVertex, newNormals] = SmoothLaplace(v, p, normals, iterations, weight, noderestrictions, max_displacement)

c=p(:,1);
a=p(:,2:end);
ind=(a==0);
a1=a+ind;
ind2=(a~=0);
b=reshape(a1,size(a1,1),size(a1,2));
noderestrictions_indexe=find(noderestrictions==1);
org=v;
for i=1:iterations  
  x=v(b,1); x=reshape(x,size(b,1),size(b,2));
  y=v(b,2); y=reshape(y,size(b,1),size(b,2));
  z=v(b,3); z=reshape(z,size(b,1),size(b,2));
  
  x=x-x.*ind; x=sparse(x);
  y=y-y.*ind; y=sparse(y);
  z=z-z.*ind; z=sparse(z);
  
  x1=ind2.*kron(ones(1,size(x,2)),v(:,1));
  y1=ind2.*kron(ones(1,size(y,2)),v(:,2));
  z1=ind2.*kron(ones(1,size(z,2)),v(:,3));
  
  d1=sum(x-x1,2)./c;
  d2=sum(y-y1,2)./c;
  d3=sum(z-z1,2)./c;
  displacement=[d1 d2 d3];
  
  x=normals(b,1); x=reshape(x,size(b,1),size(b,2));
  y=normals(b,2); y=reshape(y,size(b,1),size(b,2));
  z=normals(b,3); z=reshape(z,size(b,1),size(b,2));
  
  x=x-x.*ind; x=sparse(x);
  y=y-y.*ind; y=sparse(y);
  z=z-z.*ind; z=sparse(z);
  
  x1=ind2.*kron(ones(1,size(x,2)),normals(:,1));
  y1=ind2.*kron(ones(1,size(y,2)),normals(:,2));
  z1=ind2.*kron(ones(1,size(z,2)),normals(:,3));
  
  d1=sum(x-x1,2)./c;
  d2=sum(y-y1,2)./c;
  d3=sum(z-z1,2)./c;
  
  displacement_normals=[d1 d2 d3];
  tmp=v+displacement*weight;
  tmp=tmp-org; tmp=sum(tmp.^2,2).^(0.5);
  tmp=find(tmp>max_displacement);
  if ~isempty(tmp)
     noderestrictions(tmp)=1; 
     noderestrictions_indexe=find(noderestrictions==1);
  end
  displacement(noderestrictions_indexe,:)=0;
  displacement_normals(noderestrictions_indexe,:)=0; 
  v=v+displacement*weight;
  normals=normals+displacement_normals*weight;
  disp(['Triangle surface smoothing step ' num2str(i) ' finished!']);
end
newVertex=v;
newNormals=normals;
