clear variables
close all
home

filename = 'test.conductivity.stoch.h5';
%file_id = H5F.open(filename);
h5disp(filename);
data = h5read(filename,'/stochastic');

%order of data is x/y/z

data(:,:,[1,2,9,10]) = -6;
data(:,:,3:8) = -4;

h5write(filename,'/stochastic',data);

