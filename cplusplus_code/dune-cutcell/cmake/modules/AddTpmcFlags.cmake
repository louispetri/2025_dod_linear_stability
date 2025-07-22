#
# Module providing convenience methods for compile binaries with Tpmc support.
#
# Provides the following functions:
#
# add_dune_tpmc_flags(target1 target2 ...)
#
# adds Tpmc flags to the targets for compilation and linking
#
function(add_dune_tpmc_flags)
  if(TPMC_FOUND)
    include(CMakeParseArguments)
    cmake_parse_arguments(ADD_TPMC "OBJECT;SOURCE_ONLY" "" "" ${ARGN})
    if(ADD_TPMC_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      include_directories(${TPMC_INCLUDE_DIRS})
    else(ADD_TPMC_SOURCE_ONLY)
      if(NOT ADD_TPMC_OBJECT)
        foreach(_target ${ADD_TPMC_UNPARSED_ARGUMENTS})
          target_link_libraries(${_target} ${TPMC_LIBRARIES})
        endforeach(_target ${ADD_TPMC_UNPARSED_ARGUMENTS})
      endif(NOT ADD_TPMC_OBJECT)
      set(_prefix TARGET)
      set_property(${_prefix}  ${ADD_TPMC_UNPARSED_ARGUMENTS} APPEND
        PROPERTY
        COMPILE_DEFINITIONS ENABLE_TPMC)
      include_directories(${TPMC_INCLUDE_DIRS})
    endif(ADD_TPMC_SOURCE_ONLY)
  endif(TPMC_FOUND)
endfunction(add_dune_tpmc_flags)
