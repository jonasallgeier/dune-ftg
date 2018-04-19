# Generation of random fields using DUNEs random field generator and Matlab

To generate all necessary fields for your dune model, you have to do the following things
1. define __common parameters__ and properties of all fields via the `set_common_traits.m` function
2. define all field __specific parameters__ in respective `FIELDNAME.m` files
3. adapt `main.m` to include all your fields
4. run `main.m`

In general, randomfields of dune can be aquired using the following syntax:
`\path\to\fieldgenerator .\modelling.ini -output.dune \path\to\output`
modifications can be made, for example by using matlabs hdf5 high level functions.

