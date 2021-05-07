function [Arrows, Balls] = coilspace_geo2sci(path,metric,parameter,scale)
scaling=1;
if(nargin>3)
   scaling=scale; 
end

f=fopen(path,'r');
line=fgetl(f);
line=fscanf(f,'%s\n');
fclose(f);
line=line(1:end-2);
a=strsplit(line,';');

n=zeros(length(a),3);
f=zeros(length(a),2);
d=zeros(length(a),1);

count=0;
count2=1;
r=zeros(length(a),3);
m=zeros(length(a),1);
for i=1:length(a)
    tmp=a{i};
    if (~isempty(tmp)) 
     if strcmp(tmp(1:2),'VP')
      u=textscan(tmp, 'VP(%f,%f,%f){%f,%f,%f};');
      count=count+1;
      n(count,:)=cell2mat(u(1:3));
      v=cell2mat(u(4:6));
      mag=norm(v);
      v=v./mag;
     
      go=1;
      if (nargin>2)
       if (strcmp(parameter,'mirrored'))
         go=[1 -1];  
       end
      end
      ind=size(n,1);
      for j=1:length(go)
       s=[count count+1];
       f(count,:)=s;
       d(count,:)=mag;
       
       p=n(count,:)+go(j)*v*scaling;
       count=count+1;
       n(count,:)=p;
      end
     end
     if strcmp(tmp(1:2),'SP')
      u=textscan(tmp, 'SP(%f,%f,%f){%f};');
      r(count2,:)=cell2mat(u(1:3));
      m(count2)=cell2mat(u(4));
      count2=count2+1;
     end
    end
end
n=n(1:count,:);
f=f(1:2:count-1,:);
m=m(1:count2-1);
r=r(1:count2-1,:);
d=d(1:2:count-1);
Arrows.node=n';
Arrows.line=f';
Arrows.field=d';

Balls=[];
if nargin>1   
 [u a b]=unique(r,'rows','stable');
 Balls.node=u';

 d=zeros(length(u),1);
 for i=1:length(u)
  ind=find(b==i);
  tmp=m(ind);
  if (strcmp(metric,'mean'))
   c=mean(tmp);
  elseif (strcmp(metric,'median'))
   c=median(tmp);  
  elseif (strcmp(metric,'max'))
   c=max(tmp);
  elseif (strcmp(metric,'min'))
   c=min(tmp);
  end
  d(i,1)=c;
 end
 Balls.field=d';
end



