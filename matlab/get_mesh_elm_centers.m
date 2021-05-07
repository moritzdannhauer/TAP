function r = get_mesh_elm_centers(Field)

n=Field.node';
c=Field.cell';
if size(c,1)>size(c,2)
    c=c';
end
d=reshape(n(c(:),:)',12,size(c,2))';
r=zeros(length(d),3);
r(:,1)=mean([d(:,1) d(:,4) d(:,7) d(:,10)],2);
r(:,2)=mean([d(:,1+1) d(:,4+1) d(:,7+1) d(:,10+1)],2);
r(:,3)=mean([d(:,1+2) d(:,4+2) d(:,7+2) d(:,10+2)],2);


