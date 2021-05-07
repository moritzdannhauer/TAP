function [Arrows, Balls] = geo2scirun(path,metric,parameter,scale)
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

n=[];
f=[];
d=[];

r=[];
m=[];
for i=1:length(a)
    tmp=a{i};
    if (~isempty(tmp)) 
     if strcmp(tmp(1:2),'VP')
      u=textscan(tmp, 'VP(%f,%f,%f){%f,%f,%f};');
      n=[n; cell2mat(u(1:3))];
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
       p=n(ind,:)+go(j)*v*scaling;
       n=[n; p];
       s=[length(n)-1 length(n)];
       f=[f; s];
       d=[d; mag];
      end
     end
     if strcmp(tmp(1:2),'SP')
      u=textscan(tmp, 'SP(%f,%f,%f){%f};');
      r=[r; cell2mat(u(1:3))];
      m=[m; u(4)];
     end
    end
end
Arrows.node=n';
Arrows.line=f';
Arrows.field=d';

Balls=[];
if nargin>1   
 [u a b]=unique(r,'rows','stable');
 Balls.node=u';

 d=[];
 for i=1:length(u)
  ind=find(b==i);
  tmp=cell2mat(m(ind));
  if (strcmp(metric,'mean'))
   c=mean(tmp);
  elseif (strcmp(metric,'median'))
   c=median(tmp);  
  elseif (strcmp(metric,'max'))
   c=max(tmp);
  elseif (strcmp(metric,'min'))
   c=min(tmp);
  end
  d=[d; c];
 end
 Balls.field=d';
end




