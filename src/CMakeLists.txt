add_executable("transient" transient.cc)
add_executable("moments" moments.cc)
target_link_dune_default_libraries("transient")
target_link_dune_default_libraries("moments")
dune_symlink_to_source_files(FILES modelling.ini boundary.groundwater boundary.geoelectrics boundary.transport boundary.moments_c0 electrode.configuration well.configuration fields results)

