clear variables
clc

name = 'conductivity';


subfields.ztop = [10 18 20]-20;
subfields.number = length(subfields.ztop);
subfields.mean = [-18.421 -6.9077 -13.816];
subfields.corrlength = [30 2 10];
subfields.model = {'exponential' 'exponential' 'exponential'};
subfields.variance = [0.5 2 1];
subfields.transform = 'logNormal';

% subfields.ztop = [10 18 20];
% subfields.number = length(subfields.ztop);
% subfields.mean = [-18.421 -6.9077 -13.816];
% subfields.corrlength = [4 10 2];
% subfields.model = {'exponential' 'exponential' 'exponential'};
% subfields.variance = [0 1 0.5];
% subfields.transform = 'logNormal';

save(name)
