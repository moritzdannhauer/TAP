function [mindis index dis] = min_distance(p, points)

index=zeros(size(points,1),1);
mindis=zeros(size(points,1),1);
dis=zeros(size(points,1),1);
bool=0;
for i=1:size(points,1)
    q=p;
    q(:,1)=q(:,1)-points(i,1);
    q(:,2)=q(:,2)-points(i,2);
    q(:,3)=q(:,3)-points(i,3);
    q=sum(q.^2,2).^(0.5);
    dis=q;
    tmp=find(min(q)==q);
    if length(tmp)>1
       bool=1; 
    end
    index(i,1)=tmp(1);
    mindis(i,1)=q(index(i,1));
end

if bool==1
    disp('at least for one point the minimal distance could not uniquely be determined ');
end
