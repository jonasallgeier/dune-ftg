add_executable("transient" transient.cc)
add_executable("moments" moments.cc)
add_executable("print" print.cc)
target_link_dune_default_libraries("transient")
target_link_dune_default_libraries("moments")
target_link_dune_default_libraries("print")
dune_symlink_to_source_files(FILES modelling.ini electrode.configuration well.configuration fields results boundary)

