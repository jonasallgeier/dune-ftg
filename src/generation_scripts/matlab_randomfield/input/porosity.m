clear variables
clc

name = 'porosity';


subfields.ztop = [20]-20;
subfields.number = length(subfields.ztop);
subfields.mean = [.25];
subfields.corrlength = [1];
subfields.model = {'exponential'};
subfields.variance = [0];
subfields.transform = 'none';

save(name)
