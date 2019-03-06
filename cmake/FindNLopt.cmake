# - Find NLopt
# Find the native NLopt includes and library
#
# NLOPT_LIBRARIES - Target for NLopt
# NLOPT_FOUND     - True if nlopt found.

# Erase what has been found in favor of local version
if(NOT TARGET nlopt)
	set(NLOPT_CXX           ON  CACHE BOOL "enable cxx routines"             FORCE)
	set(BUILD_SHARED_LIBS   OFF CACHE BOOL "Build NLopt as a shared library" FORCE)
	set(NLOPT_PYTHON        OFF CACHE BOOL "build python bindings"           FORCE)
	set(NLOPT_OCTAVE        OFF CACHE BOOL "build octave bindings"           FORCE)
	set(NLOPT_MATLAB        OFF CACHE BOOL "build matlab bindings"           FORCE)
	set(NLOPT_GUILE         OFF CACHE BOOL "build guile bindings"            FORCE)
	set(NLOPT_SWIG          OFF CACHE BOOL "use SWIG to build bindings"      FORCE)
	set(NLOPT_LINK_PYTHON   OFF CACHE BOOL "link Python libs"                FORCE)
	add_subdirectory(${THIRD_PARTY_DIR}/nlopt nlopt EXCLUDE_FROM_ALL)
	target_compile_options(nlopt PRIVATE "-Wno-deprecated-declarations")
endif()
set(NLOPT_LIBRARIES nlopt)

# handle the QUIETLY and REQUIRED arguments and set NLopt_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT DEFAULT_MSG NLOPT_LIBRARIES)
