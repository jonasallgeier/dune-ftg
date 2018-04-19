clear variables
close all
home

load common_traits.mat

names = {'conductivity','storativity','porosity'};

% definition of general properties of grid
extents = [100 70 20]; % X, Y, Z
cells = [50 35 10]; % X, Y, Z


filename = strcat(outputfolder(2:end),outputtrunk,'.fieldList');

fileID = fopen(filename,'w');

fprintf(fileID,'%s\n','[ grid ]');
fprintf(fileID,'%s "%.2f %.2f %.2f"\n','extensions =',extents);
fprintf(fileID,'%s "%d %d %d"\n','cells =',cells);
fprintf(fileID,'\n');

% print field names
fprintf(fileID,'%s\n','[ randomField ]');
thestring = strcat(names(1));
for i = 2:length(names)
    thestring = strcat(thestring," ",names(i));
end
fprintf(fileID,'%s "%s"\n','types =',thestring);
fprintf(fileID,'\n');

fclose(fileID);

conductivity
field_generation_routine
porosity
field_generation_routine
storativity
field_generation_routine
