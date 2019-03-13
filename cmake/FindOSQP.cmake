# - Find OSQP
# Find the native OSQP includes and library
#
# OSQP_LIBRARIES - Target for OSQP
# OSQP_FOUND     - True if OSQP found.

# Erase what has been found in favor of local version
if(NOT TARGET osqp)
	# set(PRINTING OFF CACHE BOOL "print OSQP information" FORCE)
	add_subdirectory(${THIRD_PARTY_DIR}/osqp osqp EXCLUDE_FROM_ALL)
	target_compile_options(osqp PRIVATE "-Wno-deprecated-declarations")
endif()
set(OSQP_LIBRARIES osqp)

# handle the QUIETLY and REQUIRED arguments and set OSQP_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OSQP DEFAULT_MSG OSQP_LIBRARIES)
