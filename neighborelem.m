function [vert_neigh, vert_tri_neigh] = neighborelem(v,f)

%matlabpool close force
%matlabpool open
nrv=size(v,1);
nrf=size(f,1);
[t,ind] = sort(f(:));
[uni uni_count]=unique_count(t);

[I,J] = ind2sub(nrf,ind);
ende=cumsum(uni_count);
begin=[1; ende+1];

vert_tri_neigh=cell(nrf,1);
vert_neigh=cell(nrf,1);
max_elem=max(uni_count(:))*2+1;
t1=zeros(max_elem,1);

%parfor i=1:nrv
for i=1:nrv
  tmp=I(begin(i):ende(i));
  t2=t1;
  t2(1:max(size(tmp))+1)=[size(tmp,1) tmp'];
  vert_tri_neigh{i}=t2;
  tmp=setdiff(unique(f(tmp,:)),uni(i));
  nr2=size(tmp,1);
  t2=t1;
  if(size(tmp,1)<size(tmp,2))
    tmp=tmp';  
  end
  t2(1:max(size(tmp))+1)=[nr2 tmp'];
  vert_neigh{i}=t2;
 % progress_bar(i,nrv,0.01);
end

vert_tri_neigh=reshape(cell2mat(vert_tri_neigh), max_elem, nrv)';
vert_neigh=reshape(cell2mat(vert_neigh),max_elem, nrv)';

vert_neigh=sparse(vert_neigh(:,1:max(vert_neigh(:,1))+1));
vert_tri_neigh=sparse(vert_tri_neigh(:,1:max(vert_tri_neigh(:,1))+1));

%matlabpool close



