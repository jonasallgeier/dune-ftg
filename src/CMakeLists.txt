add_executable("dune-ftg" dune-ftg.cc)
target_link_dune_default_libraries("dune-ftg")
dune_symlink_to_source_files(FILES modelling.ini boundary.groundwater boundary.geoelectrics boundary.transport conductivity.field storativity.field porosity.field electrode.configuration2D electrode.configuration3D test.conductivity.field test.porosity.field test.storativity.stoch.h5 test.conductivity.stoch.h5 test.porosity.stoch.h5 test.storativity.trend test.conductivity.trend test.porosity.trend test.storativity.xdmf test.conductivity.xdmf test.porosity.xdmf test.fieldList test.storativity.field)
