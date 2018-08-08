clear variables
clc

name = 'conductivity';


subfields.ztop = [10 18 20]-20;
subfields.number = length(subfields.ztop);
subfields.mean = [-18.421 -6.9077 -13.816];
%subfields.mean = [-11.512925 -11.512925 -11.512925];
subfields.corrlength = [30 2 10];
subfields.model = {'exponential' 'exponential' 'exponential'};
subfields.variance = [0.5 2 1];
%subfields.variance = [0 0 0];
subfields.transform = 'logNormal';


save(name)
