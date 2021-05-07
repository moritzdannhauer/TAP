function make_brainsight_files(subject, target, A, outputfolder, ht, REVERSE_COIL_CURRENT)
  point=A(1:3,4);
  if (~exist([subject '/' outputfolder ]))
       mkdir([subject '/' outputfolder ]);
  end
  f=fopen([subject '/' outputfolder '/' 's' subject '_' target '_hair' sprintf('%.1f',ht) 'mm.txt' ],'w');
  fprintf(f,'# Version: 12\n');
  fprintf(f,'# Coordinate system: Brainsight\n');
  fprintf(f,'# Created by: Brainsight 2.4.6\n');
  fprintf(f,'# Units: millimetres, degrees, milliseconds, and microvolts\n');
  fprintf(f,'# Encoding: UTF-8\n');
  fprintf(f,'# Notes: Each column is delimited by a tab. Each value within a column is delimited by a semicolon.\n');
  fprintf(f,'# Sample\tName\tSession\tName\tIndex\tAssoc. Target\tLoc. X\tLoc. Y\tLoc. Z\tm0n0\tm0n1\tm0n2\tm1n0\tm1n1\tm1n2\tm2n0\tm2n1\tm2n2\tDist. to Target\tTarget Error\tAngular Error\tTwist Error\tDate\tTime\tCreation Cause\tCrosshairs Driver\tOffset\tComment\n');
  datum=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
  datum_1=datum(1:strfind(datum,' ')-1);
  datum_2=datum(strfind(datum,' ')+1:end);
  %if (REVERSE_COIL_CURRENT==0)
    fprintf(f,['Subject' '\t1\tSession\t2\t1\t' subject target 'HT' num2str(ht) 'mm' ' \t' num2str(point(1)) '\t'  num2str(point(2)) '\t' num2str(point(3)) '\t' num2str(A(1,1),4) '\t' num2str(A(2,1),4) '\t' num2str(A(3,1),4) '\t'  num2str(A(1,2),4) '\t' num2str(A(2,2),4) '\t' num2str(A(3,2),4) '\t' num2str(A(1,3),4) '\t' num2str(A(2,3),4) '\t' num2str(A(3,3),4) '\t' '0.0000\t0.0000\t0.0000\t0.0000\t' datum_1 '\t' datum_2 '.097\tButton\tMouse\t0.0000\t(null)\n']);
  %else
  %  fprintf(f,['Subject' '\t1\tSession\t2\t1\t' subject  target 'HT' num2str(ht) 'mm' '_REVERSECOILCURRENT' ' \t' num2str(point(1)) '\t'  num2str(point(2)) '\t' num2str(point(3)) '\t' num2str(A(1,1),4) '\t' num2str(A(2,1),4) '\t' num2str(A(3,1),4) '\t'  num2str(A(1,2),4) '\t' num2str(A(2,2),4) '\t' num2str(A(3,2),4) '\t' num2str(A(1,3),4) '\t' num2str(A(2,3),4) '\t' num2str(A(3,3),4) '\t' '0.0000\t0.0000\t0.0000\t0.0000\t' datum_1 '\t' datum_2 '.097\tButton\tMouse\t0.0000\t(null)\n']);  
  %end
  fclose(f);
 
disp(['wrote file: ' outputfolder '/s' subject '_' target '_hair' sprintf('%.1f',ht) 'mm.txt' ]);
