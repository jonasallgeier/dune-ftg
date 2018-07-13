# dune-ftg
## a model for groundwater flow, solute transport and ERT/geoelectrics in the subsurface

_This project is currently in development, as indicated by the version number (0.#.#). Do not expect this code to fully work before a release version (1.0.0)._

This code aims to simulate stationary groundwater flow in a confined aquifer, advective-dispersive solute transport depending on this groundwater flow and finally a stack of geoelectrics models simulating electrical resistivity tomography (ERT) of the subsurface. The ERT models are able to convert concentration data of the transport model to electrical conductivity in order to simulate salt tracer tests. A linear relationship is used for the transformation (Archies law). In addition to a transient formulation, moment-generating PDEs, as described by <cite> Pollock & Cirpka (2010)</cite>, are used to directly evaluate the temporal moments of solute concentration and electrical potential perturbation.

The implementation is based on a parallelized cell-centered finite volume method used on rectilinear three-dimensional grids, which do not have to be equidistant.

This is a DUNE model, meaning that this portion of code depends heavily on several other modules of the DUNE project. Information regarding DUNE can be obtained from [the dune project homepage][dune]. This module depends directly on `dune-modelling` (tested with version 1.0) and `dune-randomfield` (tested with version 1.0), who have their own dependencies. One of these is dune-pdelab (tested with version 2.5-dev).

The module can be compiled as any other DUNE module. It comprises three different executables, which can be used via the command line after compilation.


* `src/print`: save parameter fields in subsampled vtk format for visualization
* `src/transient`: a transient model for flow, transport and ERT during salt tracer tests
* `src/moments`: a moment-generating model for flow transport and ERT during salt tracer tests

User specified files include:

* `src/boundary/*`: boundary condition files; one for each model type
* `src/*.configuration`: files for the setup of electrodes and groundwater wells
* `src/fields/*`: a set of hdf5 property fields in `dune-randomfield` format
* `src/modelling.ini`: a text file containing all run time relevant model information

Horizontally layered property fields can be generated using the provided matlab wrapper for `dune-randomfield`, which is located in `src/generation_scripts/matlab_randomfield`.

The actual source code is located in `dune/ftg/` and `dune/ftg/override`. While the former contains all model-specific classes, the latter comprises modified classes of `dune-pdelab`, `dune-randomfield` and `dune-modelling`. It might be possible that future versions of those modules might be in conflict with these overloaded files. In doubt, check the source code. This is encouraged anyway.

For further information regarding this specific dune module, feel free to contact me.

--------------------------------------------------------------------------

Pollock, D., & Cirpka, O. A. (2010). Fully coupled hydrogeophysical inversion of synthetic salt tracer experiments. Water Resources Research, 46(7).

[dune]: https://dune-project.org/


