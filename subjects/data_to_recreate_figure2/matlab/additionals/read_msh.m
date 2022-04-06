function [tri_mesh, tet_mesh, tri_mesh_ef, tet_mesh_ef] = read_msh(filename)
%this only works for ascii only, for tet and tri elements; as well as gmsh *.pos files
tags=[2 4];
nelem=[3 4];

f=fopen(filename,'r');
line = fgetl(f);
line = fgetl(f);
line = fgetl(f);

nodes_present=0;
elem_present=0;
data_present=0;
tri_mesh_ef=[]; 
tet_mesh_ef=[];
tri_mesh=[]; 
tet_mesh=[];

while ~feof(f)
line = fgetl(f);  
if (strcmp(line(1),'$') && ~strcmp(line(2:4),'End'))
 if (strcmp(line,'$Nodes'))
  nn=sscanf(fgetl(f),'%d\n');
  if(nn<1)
   disp('Something is wrong! Number of mesh nodes is smaller than 1?');
   fclose(f);
   tri_mesh=[];tet_mesh=[];
   return;
  end
  nodes=zeros(nn,3);
  for i=1:nn
   line = sscanf(fgetl(f),'%d %f %f %f\n');
   nodes(i,:)=line(2:4)';
  end
  line = fgetl(f);
  if (strcmp(line,'$EndNodes'))
     nodes_present=1;  
  else
     disp('Number nodes does not match up (end $EndNodes missing)');
     fclose(f);
     tri_mesh=[];tet_mesh=[];
     return;
  end
 else
 if (strcmp(line,'$Elements')) %Elements?
  ne=sscanf(fgetl(f),'%d\n');
  if(ne<1)
   disp('Something is wrong! Number of mesh elements is smaller than 1?');
   fclose(f);
   tri_mesh=[];tet_mesh=[];
   return;
  end
  elements=zeros(ne,7); % 7= max. 4 element nodes + 1 label + 1 element type + element number
  for i=1:ne
   line = str2num(fgetl(f));
   elements(i,1) = line(1);
   elements(i,2) = line(2);
   ind=3 + line(3);
   elements(i,3) = line(ind); %last of those should be the label
   tmp=nelem(find(elements(i,2)==tags))-1;
   elements(i,4:4+tmp) = line(ind+1:ind+1+tmp);
  end
  line = fgetl(f);
  if (strcmp(line,'$EndElements'))
    elem_present=1;  
  else 
    disp('Number elements does not match up (end $EndElements missing)'); 
    fclose(f);
    tri_mesh=[];tet_mesh=[];
    return;
  end
else
 if(strcmp(line,'$ElementData'))
  line = fgetl(f); %read in headers
  line = fgetl(f);
  line = fgetl(f);
  line = fgetl(f);
  line = fgetl(f);
  line = fgetl(f);
  line = fgetl(f);
  line = fgetl(f); % components? assuming a maximum of 3 components
  nr_comp=str2num(line);
  line = fgetl(f); % number of data values
  nr_lines=str2num(line);
  if(nr_comp>=1 && nr_lines>0) 
    data=zeros(nr_lines,nr_comp+1);
    if (nr_comp==1)
       for i=1:nr_lines
          line=sscanf(fgetl(f),'%d %f\n'); 
          data(i,1:2)=line(1:2);
       end
    else
       if (nr_comp==2)
        for i=1:nr_lines   
         line=sscanf(fgetl(f),'%d %f %f\n'); 
         data(i,1:nr_comp+1)=line(1:3); 
        end
       else
         if (nr_comp==3)
          for i=1:nr_lines   
           line=sscanf(fgetl(f),'%d %f %f %f\n'); 
           data(i,1:nr_comp+1)=line(1:4); 
          end
         else
          if(nr_comp==9)
           for i=1:nr_lines   
            line=sscanf(fgetl(f),'%d %f %f %f %f %f %f %f %f %f\n'); 
            data(i,1:nr_comp+1)=line(1:10); 
           end  
          else    
           disp('Number of data vector component can be 1-3 or 9.'); 
           fclose(f);
           tri_mesh=[];tet_mesh=[];
           return;
          end
         end
       end
    end
   end
   line = fgetl(f);
   if (strcmp(line,'$EndElementData'))
    data_present=1;
   else
     disp('The definition of element data is inconsistent.');   
     fclose(f);
     tri_mesh=[];tet_mesh=[];
     return;
   end
   end   
   end  
  end 
end

%split elements to outputs
if ((nodes_present && elem_present && data_present) || feof(f))
   if (exist('data','var'))
    [~, resort]=sort(data(:,1));
    data=data(resort,2:end);
   end
   [~, resort]=sort(elements(:,1));
   elements=elements(resort,2:end);
   ind_tri=find(elements(:,1)==tags(1)); %triangles
   ind_tet=find(elements(:,1)==tags(2)); %tetrahedra
   tet_labels=elements(ind_tet,2); 
   tri_labels=elements(ind_tri,2);
   tri_elem=elements(ind_tri,3:3+nelem(1)-1);
   tet_elem=elements(ind_tet,3:3+nelem(2)-1);   
   if (min(size(tri_elem))~=0)
    tri_mesh_ef.node=nodes';
    tri_mesh_ef.face=tri_elem';
    if (data_present==1 && length(ind_tri)>1 && length(tet_elem)+length(tri_elem)==length(data))
      tri_mesh_ef.field=data(ind_tri,:)';  
    end
    tri_mesh.node=tri_mesh_ef.node;
    tri_mesh.face=tri_mesh_ef.face;
    tri_mesh.field=tri_labels';
   else
    tri_mesh_ef=[]; 
   end
   if (min(size(tet_elem))~=0)
    tet_mesh_ef.node=nodes';
    tet_mesh_ef.cell=tet_elem';
    if (data_present==1 && length(ind_tet)>1 && length(tet_elem)+length(tri_elem)==length(data))
      tet_mesh_ef.field=data(ind_tet,:)';
    end
     tet_mesh.node=tet_mesh_ef.node;
     tet_mesh.cell=tet_mesh_ef.cell;
     tet_mesh.field=tet_labels';
   else
    tet_mesh=[];
   end
   return;
end 
  
end

if (~(nodes_present && elem_present))
   fclose(f);
   tri_mesh_ef=[];tet_mesh=[];
   return;   
end




