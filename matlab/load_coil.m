function [coil, didt] = load_coil(path)
maxi=1e6; %a coil discretized with 1e6 dipoles or file lines seems an upper bound
end_count=0;
count=1;
f=fopen(path,'r');
nodes=[];
didt=nan;
disp(['[TAP] Read in ' path ' file. ']);
while ~feof(f)
  end_count=end_count+1;  
  line=fgetl(f);
  if (line(1)~='#')
    a=strsplit(line,' '); 
    if (length(a)==1)
      nodes=zeros(str2num(cell2mat(a)),3); 
      vec=zeros(str2num(cell2mat(a)),3); 
    else
      if (~isempty(nodes))
          if length(nodes)>1 && length(a)==6
            nodes(count,1)=str2num(a{1});
            nodes(count,2)=str2num(a{2});
            nodes(count,3)=str2num(a{3});
              vec(count,1)=str2num(a{4});
              vec(count,2)=str2num(a{5});
              vec(count,3)=str2num(a{6});
            count=count+1;
          end
      end
    end
  else
   if count==1
    tmp=char(extractBetween(line,'dIdtmax=',';')); 
    if (~isempty(tmp))
      didt=str2num(tmp);
    end
   end
  end
  if (end_count>=maxi)
     disp('Internal error: Coil was not properly read in. Coil transformation matrix generation questionable!'); 
     fclose(f);
     break;
  end
end
fclose(f);
coil.node=nodes';
coil.field=vec';
