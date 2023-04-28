function [optCoil] = bestcoil_geo2sci(path)

f=fopen(path,'r');
line=fgetl(f);
line=fscanf(f,'%s\n');
fclose(f);
line=line(1:end-2);
a=strsplit(line,';');
r=zeros(length(a),3);
m=zeros(length(a),1);
count=0;
for i=1:length(a)
    tmp=a{i};
    if (~isempty(tmp)) 
     if strcmp(tmp(1:2),'SP')
      u=textscan(tmp, 'SP(%f,%f,%f){%f};');
      r(i,:)=cell2mat(u(1:3));
      m(i)=cell2mat(u(4));
      count=count+1;
     end
    end      
end
r=r(1:count,:);
m=m(1:count);
optCoil.node=r';
optCoil.field=m'; %m is the mag. dipole vector only


