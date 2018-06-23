clear variables
close all
home

cm.number_of_processors = 2;
field_resolution = 'intermediate';

% directory information
cm.outputfolder = ' ./'; % space in front is necessary!
cm.outputtrunk = 'output/intermediate'; % all field files have this trunk in front of their name; it is not recommended to leave this empty
cm.input = './input'; % all field files have this trunk in front of their name; it is not recommended to leave this empty

cm.fieldgenerator_location = '../../../../../release-build/dune-randomfield/src/fieldgenerator';
cm.vtk_out = true;

cm.isequidistant = false;
cm.seed = 42; % this is a matlab rng seed, which will later deliver different pseudo random reproducible seeds for dune-randomfield
rng(cm.seed);

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
     switch field_resolution
        case 'intermediate'
            cm.vector_x = [-50:4:-34 -32 -30:-12 -11.5:.5:11.5 12:30 32 34:4:50];
            cm.vector_y = [-35:4:-19 -17 -15:-4 -3:.5:3 4:15 17 19:4:35];
            cm.vector_z = [-20:2:-12 -11:1:0];
        case 'production'
            cm.vector_x=[-50 -45 -40 -37.5 -35 -32.5 [-30:-20] [-19.5:0.5:-8] [-7.75:0.25:-4] [-3.9:.1:3.9] [4:.25:7.75] [8:0.5:19.5] [20:30] 32.5 35 37.5 40 45 50];
            cm.vector_y=[-35 -31 -27 -24 -22 [-20:-10] [-9.5:.5:-5] [-4.75:0.25:-0.5] [-.4:.1:.4] [0.5:0.25:4.75] [5:.5:9.5] [10:20] 22 24 27 31 35]; 
            cm.vector_z=[-20 -15 -12.5 -11.5 -10.7 -10.3 [-10:.1:-2] [-1.8:.2:0]];
        case 'cirpka'
            cm.vector_x=[-50 -45 -40 -37.5 -35 -32.5 [-30:-20] [-19.5:0.5:-8] [-7.7:0.3:-5] [-14:15]/3 [5.3:.3:8] [8.5:.5:15] [16:25] [27.5:2.5:40] 45 50];
            cm.vector_y=[-35 -30 -25 -22.5 [-20:-10] [-9.5:.5:-5] [-14:-1]/3 [-.2:.1:.2] [1:15]/3 [5.5:.5:10] [11:20] 22.5 25 30 35];
            cm.vector_z=[ -20 -15 -13 -11.5 -10.7 -10.3 -10.1 [-10:.1:-2] [-1.8:.2:0]];
     end
    cm.extents = [range(cm.vector_x), range(cm.vector_y), range(cm.vector_z)];
    cm.cells = [length(cm.vector_x)-1, length(cm.vector_y)-1, length(cm.vector_z)-1];
end


cm.names = {'conductivity','sigma_bg'};

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