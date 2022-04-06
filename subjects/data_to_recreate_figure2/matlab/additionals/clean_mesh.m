function Field1 = clean_mesh(n,c,f,ind)

if size(c,1)>size(c,2)
    c=c';
end

if size(n,1)<size(n,2)
   n=n';    
end

%f=ones(max(size(c)),1);

ind=setdiff(1:max(size(c)),ind);
f=f(ind);
c=c(:,ind);

%remove same elements
d=reshape(n(c(:),:)',12,size(c,2));
d=d(:);
d=reshape(d,12,size(c,2))';
[e t j]=unique(d,'rows');

%generate unique nodes
p=reshape(e',3,size(e,1)*size(e,2)/3)';
%n=size(e,1)*size(e,2)/3;
[Y,I,J] = unique(p,'rows');   %p=Y(J,:);  Y=p(I,:);
[k m]=sort(J); %k=J(m);

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

fin_cell(m(1))=1;
fin_cell=reshape(fin_cell,4,size(e,1))';

%label ?
f=f(t);

clear Field1;
Field1.node=Y';
Field1.field=f';
Field1.cell=uint32(fin_cell)';

