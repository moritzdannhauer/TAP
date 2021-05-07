function BS_Trans_mat = read_brainsight_file(file)

f=fopen(file,'r');
count=1;
max_count=1e6;
tag='# Sample';
while true
 line = fgetl(f);
 if(strcmp(line(1:length(tag)),tag)) 
   break;
 end
 count=count+1;
 if (count>max_count)
   break;
 end
end
BS_Trans_mat={};
count=1;
while ~feof(f)
 line = fgetl(f);
 line = regexprep(line, ' ', '\t');
 str=strsplit(line,'\t');
 if (strcmp(line(1),'#'))
     break;
 end
 Sample=str{1};
 SampleNr=str{2};
 Session=str{3};
 SessionNr=str{4};
 Index=str{5};
 Target=str{6};
 X=str2num(cell2mat(str(7)));
 Y=str2num(cell2mat(str(8)));
 Z=str2num(cell2mat(str(9)));
 m0n0=str2num(cell2mat(str(10)));
 m0n1=str2num(cell2mat(str(11)));
 m0n2=str2num(cell2mat(str(12)));
 m1n0=str2num(cell2mat(str(13)));
 m1n1=str2num(cell2mat(str(14)));
 m1n2=str2num(cell2mat(str(15)));
 m2n0=str2num(cell2mat(str(16)));
 m2n1=str2num(cell2mat(str(17)));
 m2n2=str2num(cell2mat(str(18)));
 T=zeros(4,4);
 T(1,1)=m0n0; T(2,1)=m0n1; T(3,1)=m0n2;
 T(1,2)=m1n0; T(2,2)=m1n1; T(3,2)=m1n2;
 T(1,3)=m2n0; T(2,3)=m2n1; T(3,3)=m2n2;
 T(1,4)=X;    T(2,4)=Y;    T(3,4)=Z;  T(4,4)=1;
 BS_Trans_mat{count}=T;
 count=count+1;
 if (count>max_count)
   break;
 end
end

disp('done.');
