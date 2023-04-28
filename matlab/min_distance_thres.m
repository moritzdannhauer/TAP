function [index] = min_distance_thres(p, points, thres)
%manhattan distance
if(length(p)<length(points))
  A=p;
  B=points;
else
  A=points;
  B=p;
end

ind=zeros(length(points),1);
count=1;
tic;
for i=1:length(A)
   tmp=B;
   tmp(:,1)=abs(tmp(:,1)-A(i,1));
   tmp(:,2)=abs(tmp(:,2)-A(i,2));
   tmp(:,3)=abs(tmp(:,3)-A(i,3));
   index1=find(tmp(:,1)<=thres);
   index2=find(tmp(:,2)<=thres);
   index3=find(tmp(:,3)<=thres);
   index=intersect(index1,index2);
   tmp=intersect(index,index3);
   if (~isempty(tmp))
    ind(count:count+length(tmp)-1)=tmp;
    last=length(tmp);
    count=count+last;
   end
end
toc
index=unique(ind(1:count-last));
disp('ok');


