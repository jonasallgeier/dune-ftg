clear variables
close all
home

% directory information
cm.outputfolder = ' ./'; % space in front is necessary!
cm.outputtrunk = 'output/testing_e'; % all field files have this trunk in front of their name; it is not recommended to leave this empty
cm.input = './input'; % all field files have this trunk in front of their name; it is not recommended to leave this empty

cm.fieldgenerator_location = '../../../../../release-build/dune-randomfield/src/fieldgenerator';
cm.vtk_out = true;

cm.isequidistant = true;
cm.seed = 2;

% definition of general properties of grid
if cm.isequidistant
    cm.extents = [100 70 20]; % X, Y, Z
    cm.cells = [50 35 10]; % X, Y, Z
    cm.cell_tops = linspace(0,cm.extents(3),cm.cells(3)+1);
    cm.cell_tops = cm.cell_tops(2:end);
    cm.vector_x = linspace(0,cm.extents(1),cm.cells(1)+1);
    cm.vector_y = linspace(0,cm.extents(2),cm.cells(2)+1);
    cm.vector_z = linspace(0,cm.extents(3),cm.cells(3)+1);
else
    cm.vector_x = [0 4 8 12 16 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 88 92 96 100]-50;
    cm.vector_y = [0 4 8 12 16 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 54 58 62 66 70]-35;
    cm.vector_z = [0 10 12 14 16 18 20]-20;
    %cm.vector_x = [0:10:50 52:3:61 70:10:100];
    %cm.vector_y = [0:10:30 32:7:42 50:10:70];
    %cm.vector_z = [0 10 18 20];
    cm.extents = [range(cm.vector_x), range(cm.vector_y), range(cm.vector_z)];
    cm.cells = [length(cm.vector_x)-1, length(cm.vector_y)-1, length(cm.vector_z)-1];
end


cm.names = {'conductivity','storativity','porosity','sigma_bg','kappa'};

filename = strcat(cm.outputfolder(2:end),cm.outputtrunk,'.fieldList');

fileID = fopen(filename,'w');

fprintf(fileID,'%s\n','[ grid ]');
fprintf(fileID,'%s "%.2f %.2f %.2f"\n','extensions =',cm.extents);
fprintf(fileID,'%s "%d %d %d"\n','cells =',cm.cells);
fprintf(fileID,'\n');

% print field names
fprintf(fileID,'%s\n','[ randomField ]');
thestring = strcat(cm.names{1});
for i = 2:length(cm.names)
    thestring = strcat(thestring," ",cm.names(i));
end
fprintf(fileID,'%s "%s"\n','types =',thestring);
fprintf(fileID,'\n');

fclose(fileID);

for i = 1:length(cm.names)
    field_generation_routine(cm.names{i},cm)
end