clear variables
close all
home

%limits = [-50 50; -35 35; -20 0];

% define vectors of electrode placement
x_values = 10/3*[-1.5 -.5 .5 1.5];
y_values = 10/3*[-1.5 -.5 .5 1.5];
z_values = linspace(-20,0,9);

% transformation of coordinates, here movement of origin
%x_values = x_values+50;
%y_values = y_values+35;
%z_values = z_values+20;

% get coordinate vector of electrodes
coordinates = zeros(length(x_values)*length(y_values)*length(z_values),3);
counter = 1;
for i = 1:length(x_values)
    for j = 1:length(y_values)
        for k = 1:length(z_values)
            coordinates(counter,:) = [x_values(i), y_values(j),z_values(k)];
            counter = counter+1;
        end
    end
end

% electrode format
identifier = 1:size(coordinates,1);
identifier = identifier';
surface = zeros(size(coordinates,1),1);
surface(coordinates(:,3)==0) = 1;
surface = logical(surface);

% assemble electrode file matrix EFM
EFM = [identifier, coordinates, surface];
fileID = fopen('electrode.configuration','w');
fprintf(fileID,'%d\n',length(x_values)*length(y_values)*length(z_values));
formatSpec = '%d %8.3f %8.3f %8.3f %d\n';
fprintf(fileID,formatSpec,[identifier,coordinates,surface]');
fclose(fileID);

