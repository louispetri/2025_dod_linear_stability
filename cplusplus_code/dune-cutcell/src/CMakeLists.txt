add_executable("simulations-scalar-transport" simulationsscalartransport.cc)
target_link_dune_default_libraries("simulations-scalar-transport")

configure_file(scalar-transport.ini ${CMAKE_BINARY_DIR}/scalar-transport.ini)
