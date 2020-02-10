include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "")
endif()

function(fixing_collisions_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FIXING_COLLISIONS_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FIXING_COLLISIONS_EXTERNAL}/.cache/${name}
        QUIET
        ${FIXING_COLLISIONS_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

function(fixing_collisions_download_catch2)
     fixing_collisions_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.5.0
    )
endfunction()

function(fixing_collisions_download_nlopt)
     fixing_collisions_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        463abb855d166503a4f2890379647a54a26ca678
    )
endfunction()

function(fixing_collisions_download_osqp)
     fixing_collisions_download_project(osqp
        GIT_REPOSITORY https://github.com/oxfordcontrol/osqp.git
        GIT_TAG        8949e678122d6949139644d2a95985765527535f
    )
endfunction()

function(fixing_collisions_download_fmt)
     fixing_collisions_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        5.3.0
    )
endfunction()

function(fixing_collisions_download_spdlog)
     fixing_collisions_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.3.1
    )
endfunction()

function(fixing_collisions_download_libigl)
     fixing_collisions_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        aea868bd1fc64f71afecd2c51e51507a99d8e3e5
    )
endfunction()

function(fixing_collisions_download_json)
     fixing_collisions_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json
        GIT_TAG v3.7.0
    )
endfunction()

function(fixing_collisions_download_cli11)
     fixing_collisions_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        8ecce8fd2c49f64c80e5757cb12d2fd1fa62f242
    )
endfunction()

function(fixing_collisions_download_finite_diff)
     fixing_collisions_download_project(finite-diff
        GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
        GIT_TAG        01e16fb39ecf6fe6367b932448028817becfb633
    )
endfunction()

function(fixing_collisions_download_tbb)
   fixing_collisions_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        20357d83871e4cb93b2c724fe0c337cd999fd14f
  )
endfunction()

# Etienne Vouga's CTCD Library
function(fixing_collisions_download_evctcd)
  fixing_collisions_download_project(EVCTCD
    GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
    GIT_TAG        e5fe5c9767207df5047e375fb20180a665ae186f
  )
endfunction()
