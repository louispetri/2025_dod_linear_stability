# Module that checks whether tpmc is available
#
# Accepts the following input variable
# TPMC_PREFIX: Prefix under which tpmc is installed
#
# The following variable will be set:
# TPMC_FOUND: whether tpmc is available
# TPMC_INCLUDE_DIRS: Include directories for tpmc
# TPMC_LIBRARIES: Full path to libraries needed to link
#   to tpmc
#
set(TPMC_FOUND TPMC_FOUND-NOTFOUND)

# find header in user supplied directory
find_path(TPMC_INCLUDE_DIR tpmc/marchinglut.hh
  PATHS ${TPMC_PREFIX}
  PATH_SUFFIXES include include/python tpmc/include tpmc/tpmc/include
  NO_DEFAULT_PATH
  DOC "Include directory with tpmc header files")
find_path(TPMC_INCLUDE_DIR tpmc/marchinglut.hh
  PATH_SUFFIXES include include/python tpmc/include tpmc/tpmc/include
  DOC "Include directory with tpmc header files")

# check header usability
include(CMakePushCheckState)
# store current check state
cmake_push_check_state()

# check if include file is working
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES ${TPMC_INCLUDE_DIR})
check_include_file_cxx(tpmc/marchinglut.hh TPMC_HEADER_USABLE)

# find library in custom directory
find_library(TPMC_LIBRARY
  NAMES tpmc_tables
  PATHS ${TPMC_PREFIX}
  PATH_SUFFIXES lib lib/python/tpmc/lib tpmc/lib
  NO_DEFAULT_PATH
  DOC "Full path to tpmc library.")
# find lib in default directory
find_library(TPMC_LIBRARY
  NAMES tpmc_tables
  PATH_SUFFIXES lib lib/python/tpmc/lib)

# check if library is working
include(CheckLibraryExists)
get_filename_component(TPMC_LIBRARY_PATH ${TPMC_LIBRARY} PATH)
check_library_exists(tpmc_tables table_cube2d_vertices "${TPMC_LIBRARY_PATH}" TPMC_LIB_FUNCTIONAL)

# restore check state
cmake_pop_check_state()

if(TPMC_LIB_FUNCTIONAL)
  set(TPMC_INCLUDE_DIRS ${TPMC_INCLUDE_DIR})
  set(TPMC_LIBRARIES ${TPMC_LIBRARY})
endif(TPMC_LIB_FUNCTIONAL)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "tpmc"
  DEFAULT_MSG
  TPMC_INCLUDE_DIR
  TPMC_LIBRARY
  TPMC_LIB_FUNCTIONAL
  TPMC_HEADER_USABLE
)

set(HAVE_TPMC ${TPMC_FOUND})

if (TPMC_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_TPMC=1"
    INCLUDE_DIRS "${TPMC_INCLUDE_DIRS}"
    LIBRARIES "${TPMC_LIBRARIES}")
endif(TPMC_FOUND)

mark_as_advanced(TPMC_INCLUDE_DIR TPMC_LIBRARY TPMC_LIB_FUNCTIONAL TPMC_HEADER_USABLE)
