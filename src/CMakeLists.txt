add_executable("dune-ftg" dune-ftg.cc)
target_link_dune_default_libraries("dune-ftg")
dune_symlink_to_source_files(FILES modelling.ini boundary.groundwater boundary.geoelectrics boundary.transport electrode.configuration3D fields)
