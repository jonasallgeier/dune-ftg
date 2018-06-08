function [] = field_generation_routine(name,cm)
load(strcat(cm.input,'/',name),'subfields');

%warning('hard coded z shift of 20!');
%subfields.ztop = subfields.ztop+20;

if ~cm.isequidistant
    delta_x = min(diff(cm.vector_x));
    delta_y = min(diff(cm.vector_y));
    delta_z = min(diff(cm.vector_z));
    
    x_fine = min(cm.vector_x):delta_x:max(cm.vector_x);
    y_fine = min(cm.vector_y):delta_y:max(cm.vector_y);
    z_fine = min(cm.vector_z):delta_z:max(cm.vector_z);
    
    x_fine = round(x_fine,4);
    y_fine = round(y_fine,4);
    z_fine = round(z_fine,4);
    
    if max(x_fine) ~= max(cm.vector_x) || max(y_fine) ~= max(cm.vector_y) || max(z_fine) ~= max(cm.vector_z)
        warning(['At least one extent is not an integer multiple of the smalles mesh size in the respective' ...
            ' direction. The field will still be generated, but inaccuracies are possible.']);
        
    end
    
    cm.cell_tops = z_fine(2:end);
    %subfields.finecells = NaN(size(subfields.ztop));
    %subfields.finecells(1) = find(cell_tops>=subfields.ztop(1),1);
    %for i = 2:subfields.number
    %    subfields.finecells(i) = find(cell_tops>=subfields.ztop(i),1)-sum(subfields.finecells(1:i-1));
    %end
    
    subfields.cells_in_layer = NaN(size(subfields.ztop,2),length(cm.vector_z(2:end)));
    subfields.cells_in_layer(1,:) = (cm.vector_z(2:end)<=subfields.ztop(1))';
    for i = 2:subfields.number
        subfields.cells_in_layer(i,:) = (cm.vector_z(2:end)<=subfields.ztop(i) & cm.vector_z(2:end)>subfields.ztop(i-1))';
    end
end

subfields.cells = NaN(size(subfields.ztop));
subfields.cells(1) = find(cm.cell_tops>=subfields.ztop(1),1);
for i = 2:subfields.number
    subfields.cells(i) = find(cm.cell_tops>=subfields.ztop(i),1)-sum(subfields.cells(1:i-1));
end

subfields.height = NaN(size(subfields.ztop));
for i = 2:subfields.number
    subfields.height(i) = subfields.ztop(i)-subfields.ztop(i-1);
end
if cm.isequidistant
    subfields.height(1) = subfields.ztop(1);
else
    subfields.height(1) = subfields.ztop(1)-min(cm.vector_z);   
end




for i = 1:subfields.number
    filebasename = strcat(name,num2str(i));
    filename = strcat(filebasename,'.field');
    fileID = fopen(filename,'w');
    %fileID = fopen('grid.ini','w');
    
    % print grid properties
    fprintf(fileID,'%s\n','[ grid ]');
    fprintf(fileID,'%s "%g %g %g"\n','extensions =',[cm.extents(1:2),subfields.height(i)]);
    if cm.isequidistant
        fprintf(fileID,'%s "%d %d %d"\n','cells =',[cm.cells(1:2),subfields.cells(i)]);
    else
        fprintf(fileID,'%s "%d %d %d"\n','cells =',[length(x_fine)-1,length(y_fine)-1,subfields.cells(i)]);
    end
    fprintf(fileID,'\n');
    
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
    if isempty(cm.seed)
        %command = strcat('mpirun -np ',32,num2str(cm.number_of_processors),32,cm.fieldgenerator_location,' ./',filename,' -input.seed 2 -output.dune ',cm.outputfolder,filebasename);
        command = strcat(cm.fieldgenerator_location,' ./',filename,' -input.seed 2 -output.dune ',cm.outputfolder,filebasename);
    else
        seed = floor(1e6*rand(1));
        %command = strcat('mpirun -np ',32,num2str(cm.number_of_processors),32,cm.fieldgenerator_location,' ./',filename,' -input.seed ',32,num2str(seed),' -output.dune ',cm.outputfolder,filebasename);
        command = strcat(cm.fieldgenerator_location,' ./',filename,' -input.seed ',32,num2str(seed),' -output.dune ',cm.outputfolder,filebasename);
    end
    [~,~] = system(command,'-echo');
end


parpool('local',cm.number_of_processors);

% create empty data container of size XxYxZ
field_data = NaN(cm.cells);
for i = 1:subfields.number
    filebasename = strcat(name,num2str(i));
    filename = strcat(filebasename,'.stoch.h5');
    data = h5read(filename,'/stochastic');
    data = data + subfields.mean(i);
    
    if cm.isequidistant
        if i == 1
            field_data(1:cm.cells(1),1:cm.cells(2),1:subfields.cells(i)) = data;
        else
            temp = cumsum(subfields.cells);
            field_data(1:cm.cells(1),1:cm.cells(2),1+temp(i-1):temp(i)) = data;
        end
    else
        temp = cumsum(subfields.height);
        %z_layermax = temp(i);
        
        if i== 1
            z_layermin = min(cm.vector_z);
        else
            z_layermin = temp(i-1)+min(cm.vector_z);
        end
        
        
        %mid point of cell
        %x_cell_midpoints = 0.5*delta_x+min(x_fine):delta_x:max(x_fine)-0.5*delta_x;
        %y_cell_midpoints = 0.5*delta_y+min(y_fine):delta_y:max(y_fine)-0.5*delta_y;
        z_cell_midpoints = 0.5*delta_z+min(z_fine):delta_z:max(z_fine)-0.5*delta_z;
        %left & right boundary of cell
        x_bl = x_fine(1:end-1);
        x_br = x_fine(2:end);
        y_bl = y_fine(1:end-1);
        y_br = y_fine(2:end);
        z_bl = z_fine(1:end-1);
        z_br = z_fine(2:end);
        
        
        
        x_loop_vector = cm.cells(1);
        y_loop_vector = cm.cells(2);
        z_loop_vector = find(subfields.cells_in_layer(i,:),1):1:find(subfields.cells_in_layer(i,:),1,'last'); % gets the cell indices in vertical direction
        
        help_x_max = cm.vector_x(2:end);
        help_x_min = cm.vector_x(1:end-1);
        help_y_max = cm.vector_y(2:end);
        help_y_min = cm.vector_y(1:end-1);
        help_z_max = cm.vector_z(2:end);
        help_z_min = cm.vector_z(1:end-1);
        help_cells_in_layer = length(subfields.cells_in_layer(i,:));
        
        parfor i_x = 1:x_loop_vector
            x_max = help_x_max(i_x);
            x_min = help_x_min(i_x);
            temp_matrix = NaN(length(y_loop_vector),help_cells_in_layer);
            for i_y = 1:y_loop_vector
                y_max = help_y_max(i_y);
                y_min = help_y_min(i_y);
                for i_z = z_loop_vector % i_z indicates the current cell index in vertical direction
                    z_max = help_z_max(i_z);
                    z_min = help_z_min(i_z);

                    temp_1 = (x_br-x_min)/delta_x;
                    temp_1(temp_1>1) = 1;
                    temp_1(temp_1<0) = 0;
                    
                    temp_2 = (x_max-x_bl)/delta_x;
                    temp_2(temp_2>1) = 1;
                    temp_2(temp_2<0) = 0;
                    weights_x = temp_1.*temp_2;
                    
                    temp_1 = (y_br-y_min)/delta_y;
                    temp_1(temp_1>1) = 1;
                    temp_1(temp_1<0) = 0;
                    
                    temp_2 = (y_max-y_bl)/delta_y;
                    temp_2(temp_2>1) = 1;
                    temp_2(temp_2<0) = 0;
                    weights_y = temp_1.*temp_2;
                    
                    temp_1 = (z_br-z_min)/delta_z;
                    temp_1(temp_1>1) = 1;
                    temp_1(temp_1<0) = 0;
                    
                    temp_2 = (z_max-z_bl)/delta_z;
                    temp_2(temp_2>1) = 1;
                    temp_2(temp_2<0) = 0;
                    weights_z = temp_1.*temp_2;
                    
                    weights_x = repmat(weights_x',[1,size(data,2),length(z_cell_midpoints)]);
                    weights_y = repmat(weights_y,[size(data,1),1,length(z_cell_midpoints)]);
                    weights_z = reshape(weights_z,[1,1,length(z_cell_midpoints)]);
                    weights_z = repmat(weights_z,[size(data,1),size(data,2),1]);
                    
                    weights = weights_x.*weights_y.*weights_z;
                    weights = weights(:,:,find(z_fine(2:end)>z_layermin,1):end);
                    
                    selection = data(weights>0);
                    weights = weights(weights>0);
                    
                    X_vec = reshape(selection,[1,size(selection,1)*size(selection,2)*size(selection,3)]);
                    weights_vec = reshape(weights,[1,size(weights,1)*size(weights,2)*size(weights,3)]);
                    average = sum(weights_vec.*X_vec)./sum(weights_vec);
                    
                    temp_matrix(i_y,i_z) = average;
                end
            end
            field_data(i_x,:,z_loop_vector) = temp_matrix(:,z_loop_vector);
        end
    end
    delete(strcat(filebasename,'.stoch.h5'))
    delete(strcat(filebasename,'.trend'))
    delete(strcat(filebasename,'.xdmf'))
    delete(strcat(filebasename,'.field'))
end



outputfilename = strcat(cm.outputfolder(2:end),cm.outputtrunk,'.',name,'.stoch.h5');
if exist(outputfilename, 'file') == 2
    delete(outputfilename)
end
h5create(outputfilename,'/stochastic',cm.cells);
h5write(outputfilename,'/stochastic',field_data);

outputfilename = strcat(cm.outputfolder(2:end),cm.outputtrunk,'.',name,'.field');
fileID = fopen(outputfilename,'w');
%fileID = fopen('grid.ini','w');

% print grid properties
fprintf(fileID,'%s\n','[ grid ]');
fprintf(fileID,'%s "%.2f %.2f %.2f"\n','extensions =',cm.extents);
fprintf(fileID,'%s "%d %d %d"\n','cells =',cm.cells);
fprintf(fileID,'%s "','vector_x =');
fprintf(fileID,'%g ',cm.vector_x);
fprintf(fileID,'" \n');
fprintf(fileID,'%s "','vector_y =');
fprintf(fileID,'%g ',cm.vector_y);
fprintf(fileID,'" \n');
fprintf(fileID,'%s "','vector_z =');
fprintf(fileID,'%g ',cm.vector_z);
fprintf(fileID,'" \n');
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

outputfilename = strcat(cm.outputfolder(2:end),cm.outputtrunk,'.',name,'.trend');
fileID = fopen(outputfilename,'w');
fprintf(fileID,'%s %.2f\n','mean = 0');
fclose(fileID);

if cm.vtk_out == true
    if cm.isequidistant
        vtkwrite(strcat('./',name,'.vtk'), 'structured_points', 'Conductivity', field_data)
    else
        X = NaN(size(field_data));
        Y = NaN(size(field_data));
        Z = NaN(size(field_data));
        for i_x = 1:size(field_data,1)
            for i_y = 1:size(field_data,2)
                for i_z = 1:size(field_data,3)
                    X(i_x,i_y,i_z) = 0.5*cm.vector_x(i_x)+0.5*cm.vector_x(i_x+1);
                    Y(i_x,i_y,i_z) = 0.5*cm.vector_y(i_y)+0.5*cm.vector_y(i_y+1);
                    Z(i_x,i_y,i_z) = 0.5*cm.vector_z(i_z)+0.5*cm.vector_z(i_z+1);
                end
            end
        end
        %vtkwrite(strcat('./',name,'.vtk'), 'structured_points', 'data', field_data, 'spacing', 2, 2, 2)
        vtkwrite(strcat('./',name,'.vtk'), 'structured_grid',X,Y,Z,'scalars','data', field_data)
    end
end
end
