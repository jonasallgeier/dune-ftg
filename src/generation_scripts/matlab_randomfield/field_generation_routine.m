if (sum(sum(subfields.ztop == cell_tops')) ~= subfields.number)
    number_false_ztops = subfields.number-sum(sum(subfields.ztop == cell_tops'));
    switch number_false_ztops>1
        case 1
            error(strcat('There are',32,num2str(number_false_ztops),' upper layer boundaries not aligning with the cells! Please check input.'));
        case 0
            error(strcat('There is an upper layer boundary not aligning with the cells! Please check input.'));
    end
end

subfields.cells = NaN(size(subfields.ztop));
subfields.cells(1) = find(cell_tops>=subfields.ztop(1),1);
for i = 2:subfields.number
    subfields.cells(i) = find(cell_tops>=subfields.ztop(i),1)-sum(subfields.cells(1:i-1));
end

subfields.height = NaN(size(subfields.ztop));
for i = 2:subfields.number
    subfields.height(i) = subfields.ztop(i)-subfields.ztop(i-1);
end
subfields.height(1) = subfields.ztop(1);





for i = 1:subfields.number
    filebasename = strcat(name,num2str(i));
    filename = strcat(filebasename,'.field');
    fileID = fopen(filename,'w');
    %fileID = fopen('grid.ini','w');
    
    % print grid properties TODO: print only the needed part
    fprintf(fileID,'%s\n','[ grid ]');
    fprintf(fileID,'%s "%.2f %.2f %.2f"\n','extensions =',[extents(1:2),subfields.height(i)]);
    fprintf(fileID,'%s "%d %d %d"\n','cells =',[cells(1:2),subfields.cells(i)]);
    fprintf(fileID,'\n');
    
    % print field names
    %fprintf(fileID,'%s\n','[ randomField ]');
    %thestring = strcat(name,'1');
    %for i = 2:subfields.number
    %    thestring = strcat(thestring," ",name,num2str(i));
    %end
    %fprintf(fileID,'%s "%s"\n','types =',thestring);
    %fprintf(fileID,'\n');
    
    %fclose(fileID);
    
    % print transformation type
    fprintf(fileID,'%s\n','[ randomField ]');
    fprintf(fileID,'%s %s\n','transform =',subfields.transform);
    fprintf(fileID,'\n');
    
    % print stochastic properties
    fprintf(fileID,'%s\n','[ stochastic ]');
    fprintf(fileID,'%s %.2f\n','corrLength =',subfields.corrlength(i));
    fprintf(fileID,'%s %s\n','covariance =',subfields.model{i});
    fprintf(fileID,'%s %.2f\n','variance =',subfields.variance(i));
    fprintf(fileID,'\n');
    
    % print mean properties
    fprintf(fileID,'%s\n','[ mean ]');
    fprintf(fileID,'%s %.2f\n','mean =',subfields.mean(i));
    fprintf(fileID,'%s\n','variance = 0');
    fprintf(fileID,'\n');
    fclose(fileID);
    
    % there might be compatibility issues between c++ stdlibs of matlab and the
    % system -> to make the following system command work you have to set the
    % environment variable LD_PRELOAD before starting matlab:
    % alias matlab='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21 /usr/local/bin/matlab -desktop'
    
    command = strcat(fieldgenerator_location,' ./',filename,' -output.dune ',outputfolder,filebasename);
    [status,cmdout] = system(command,'-echo');
end

% create empty data container of size XxYxZ
field_data = NaN(cells);
for i = 1:subfields.number
    filebasename = strcat(name,num2str(i));
    filename = strcat(filebasename,'.stoch.h5');
    data = h5read(filename,'/stochastic');
    data = data + subfields.mean(i);
    if i == 1
        field_data(1:cells(1),1:cells(2),1:subfields.cells(i)) = data;
    else
        temp = cumsum(subfields.cells);
        field_data(1:cells(1),1:cells(2),1+temp(i-1):temp(i)) = data;
    end
    
    delete(strcat(filebasename,'.stoch.h5'))
    delete(strcat(filebasename,'.trend'))
    delete(strcat(filebasename,'.xdmf'))
    delete(strcat(filebasename,'.field'))
end



outputfilename = strcat(outputfolder(2:end),outputtrunk,'.',name,'.stoch.h5');
h5create(outputfilename,'/stochastic',cells);
h5write(outputfilename,'/stochastic',field_data);

outputfilename = strcat(outputfolder(2:end),outputtrunk,'.',name,'.field');
fileID = fopen(outputfilename,'w');
%fileID = fopen('grid.ini','w');

% print grid properties
fprintf(fileID,'%s\n','[ grid ]');
fprintf(fileID,'%s "%.2f %.2f %.2f"\n','extensions =',extents);
fprintf(fileID,'%s "%d %d %d"\n','cells =',cells);
fprintf(fileID,'\n');

% print transformation type
fprintf(fileID,'%s\n','[ randomField ]');
fprintf(fileID,'%s %s\n','transform =',subfields.transform);
fprintf(fileID,'\n');

% print mean properties (this is needed, because the read in of randomfield requires it)
fprintf(fileID,'%s\n','[ stochastic ]');
fprintf(fileID,'%s\n','variance = 0');
fprintf(fileID,'%s\n','covariance = "exponential"');
fprintf(fileID,'\n');
fclose(fileID);

outputfilename = strcat(outputfolder(2:end),outputtrunk,'.',name,'.trend');
fileID = fopen(outputfilename,'w');
fprintf(fileID,'%s %.2f\n','mean = 0');
fclose(fileID);

if vtk_out == true
    vtkwrite(strcat('./',name,'.vtk'), 'structured_points', 'Conductivity', field_data, 'spacing', 2, 2, 2)
end