clear variables
clc

name = 'kappa';


subfields.ztop = [20]-20;
subfields.number = length(subfields.ztop);
subfields.mean = [.06];
subfields.corrlength = [1];
subfields.model = {'exponential'};
subfields.variance = [0];
subfields.transform = 'none';

save(name)