clear variables
close all
home

% directory information
outputfolder = ' ./'; % space in front is necessary!
outputtrunk = 'output/smallvariance'; % all field files have this trunk in front of their name; it is not recommended to leave this empty
fieldgenerator_location = '../../../../../release-build/dune-randomfield/src/fieldgenerator';


% definition of general properties of grid
extents = [100 70 20]; % X, Y, Z
cells = [50 35 10]; % X, Y, Z

cell_tops = linspace(0,extents(3),cells(3)+1);
cell_tops = cell_tops(2:end);

vtk_out = true;

save common_traits.mat