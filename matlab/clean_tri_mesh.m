function Field1 = clean_tri_mesh(n,c,f)

if size(c,1)>size(c,2)
    c=c';
end

if size(n,1)<size(n,2)
   n=n';    
end

d=reshape(n(c(:),:)',9,size(c,2)); %doubled triangles
d=d(:);
d=reshape(d,9,size(c,2))';
[e t j]=unique(d,'rows','stable');  

p=reshape(e',3,size(e,1)*size(e,2)/3)'; %doubled points
%n=size(e,1)*size(e,2)/3;
[Y,I,J] = unique(p,'rows','stable');
[k m]=sort(J);

fin_cell=zeros(size(k,1),1);
count=1;
for i=2:size(k,1)
  if( k(i)==k(i-1) )
      fin_cell(m(i))=count;
  else   
    count=count+1; 
    fin_cell(m(i))=count;
  end
end

fin_cell(1)=m(1);
fin_cell=reshape(fin_cell,3,size(e,1))';

f=f(t,:);

clear Field1;
Field1.node=Y';
Field1.field=f';
Field1.face=fin_cell';



