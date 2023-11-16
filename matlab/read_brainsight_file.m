function BS_Trans_mat = read_brainsight_file(file)

f=fopen(file,'r');%open file for reading
tag='# Sample';

while true %what is this while loop doing?
    line = fgetl(f);%fgetl: Read line from file, removing newline characters
    if(strcmp(line(1:length(tag)),tag)) %strcmp: string compare
        break;
    end
    if (feof(f))
        break;
    end
end
BS_Trans_mat={};
count=1;
max_count=1e6;
while ~feof(f) % feof: test for end of file; while not end of file
    line = fgetl(f);
    line = regexprep(line, 'Subject\t', 'Subject ');
    line = regexprep(line, 'Session\t', 'Session ');
    str=strsplit(line,'\t');%Split string at specified delimiter
    if (strcmp(line(1),'#')) %if the line starts w/ #
        break;
    end
    null_ind=find(strcmp(str,'(null)'));
    if (isempty(intersect(null_ind,5:16)))
        Sample=str{1};
        Session=str{2};
        Index=str{3};
        Target=str{4};
        X=str2num(cell2mat(str(5)));
        Y=str2num(cell2mat(str(6)));
        Z=str2num(cell2mat(str(7)));
        m0n0=str2num(cell2mat(str(8)));
        m0n1=str2num(cell2mat(str(9)));
        m0n2=str2num(cell2mat(str(10)));
        m1n0=str2num(cell2mat(str(11)));
        m1n1=str2num(cell2mat(str(12)));
        m1n2=str2num(cell2mat(str(13)));
        m2n0=str2num(cell2mat(str(14)));
        m2n1=str2num(cell2mat(str(15)));
        m2n2=str2num(cell2mat(str(16)));
        %making sure T elements are numeric
        check_failed = 0;
        check = {X;Y;Z;m0n0;m0n1;m0n2;m1n0;m1n1;m1n2;m2n0;m2n1;m2n2};
        for check_indx = 1:length(check)
            if isnan(check{check_indx}) || ~isnumeric(check{check_indx})
                check_failed = 1;
            end
        end
        if check_failed == 1
            disp('[TAP] One line with NaN or empty data skipped.');
        else
            T=zeros(4,4);
            T(1,1)=m0n0; T(2,1)=m0n1; T(3,1)=m0n2;
            T(1,2)=m1n0; T(2,2)=m1n1; T(3,2)=m1n2;
            T(1,3)=m2n0; T(2,3)=m2n1; T(3,3)=m2n2; %rotation of the coil
            T(1,4)=X;    T(2,4)=Y;    T(3,4)=Z;  T(4,4)=1;%center of the coil
            T(:, [1 3]) = -T(:, [1 3]); % Brainsight coordinate convention has X and Z flipped compared to the SimNIBS convention
            BS_Trans_mat{count}=T;
            count=count+1;
            if (count>max_count)
                break;
            end
        end
    end
end
disp(['[TAP] Reading ' file ' done.']);
