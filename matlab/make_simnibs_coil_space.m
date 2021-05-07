function coil = make_simnibs_coil_space(fn)

c=load(fn);

p=c(:,1:3);
f=c(:,4:6);

ind=1:3:length(p);
n=p(ind,:);
Field1.node=n';
u=[];
for i=ind
  tmp1=cross(f(i,:),f(i+1,:));
  tmp1=tmp1/norm(tmp1);
  u=[u; tmp1];
end

Field1.field=-u';

coil=Field1;

