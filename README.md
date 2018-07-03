# dune-ftg
## a model for flow, transport and geoelectrics in the subsurface

_This project is currently in development, as indicated by the version number (0.#.#). More features are planned and will be implemented. Do not expect this code to fully work before a release version (1.0.0)._

This is a DUNE model, meaning that this portion of code depends heavily on several other modules of the DUNE project. Information regarding DUNE can be obtained from [here][dune]. This module depends directly on `dune-modelling` (tested with version 1.0) and `dune-randomfield` (tested with version 1.0), who have their on dependencies. One of these is dune-pdelab (tested with version 2.5-dev).

This model comprises three different conceptual models: one for stationary groundwater flow in a confined aquifer, one for advective(-diffusive) solute transport depending on this groundwater flow and finally a stack of geoelectrics models simulating electrical resistivity tomography (ERT) of the subsurface. The ERT models are able to convert concentration data of the transport model to electrical conductivity in order to simulate salt tracer tests.

The module can be compiled as any other DUNE module. The source file is located in ./dune-ftg/src, where also all user specified data are:

* boundary condition files, one for each model type
* configuration files for electrodes and groundwater wells
* a set of hdf5 property fields in ./fields/ 
* modelling.ini, a text file containing all run time relevant model information

The compiled program can be used via the command line.

For further information regarding this specific dune-module, feel free to contact me.

[dune]: https://dune-project.org/
