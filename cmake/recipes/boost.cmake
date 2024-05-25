if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'")

# Add boost lib sources
set(BOOST_INCLUDE_LIBRARIES numeric/interval)
set(BOOST_ENABLE_CMAKE ON)

# Download and extract the boost library from GitHub
message(STATUS "Downloading and extracting boost library sources. This will take some time...")
include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE) # Needed to print downloading progress
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.7z # downloading a zip release speeds up the download
    USES_TERMINAL_DOWNLOAD TRUE 
    GIT_PROGRESS TRUE   
    DOWNLOAD_NO_EXTRACT FALSE
)
FetchContent_MakeAvailable(Boost)

# Add the boost library to the project
add_library(Boost::boost INTERFACE IMPORTED)
target_include_directories(Boost::boost INTERFACE ${boost_SOURCE_DIR})
target_compile_definitions(Boost::boost INTERFACE BOOST_ALL_NO_LIB)
target_link_libraries(Boost::boost INTERFACE Boost::numeric_interval Boost::core)
