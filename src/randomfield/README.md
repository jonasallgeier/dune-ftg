# Generation of random fields using DUNEs random field generator and Matlab

The `main.m` matlab script is able to read in horizontally defined parameter field configurations from the input directory. These are then processed using the `dune-randomfield` generator. Non-equidistant grids will create highly resolved temporary parameter fields, which are then arithmetically and volume-weighted averaged to the coarse grids.
