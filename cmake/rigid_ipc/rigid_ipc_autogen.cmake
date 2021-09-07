# Autogen helper function
function(project_autogen MAIN_TARGET PYTHON_SCRIPT TEMPLATE OUTPUT_BASE OUTPUT_SUFFIX)
    find_package(PythonInterp 3 QUIET)

    if(${PYTHONINTERP_FOUND})
    else()
        execute_process(
                COMMAND
                python -c "import sys; sys.stdout.write(sys.version)"
                OUTPUT_VARIABLE PYTHON_VERSION)

        if(NOT PYTHON_VERSION)
                MESSAGE(FATAL_ERROR "Unable to run python, ${PYTHON_SCRIPT} not running")
                return()
        endif()

        STRING(REGEX MATCH "^3\.*" IS_PYTHON3 ${PYTHON_VERSION})

        if(NOT IS_PYTHON3)
                MESSAGE(FATAL_ERROR "Unable to find python 3, ${PYTHON_SCRIPT} not running")
                return()
        else()
                SET(PYTHON_EXECUTABLE "python")
        endif()
    endif()

    execute_process(
            COMMAND
            ${PYTHON_EXECUTABLE} -c "import sympy; import sys; sys.stdout.write('ok')"
            OUTPUT_VARIABLE PYTHON_HAS_LIBS)

    if(NOT PYTHON_HAS_LIBS)
        MESSAGE(FATAL_ERROR "Unable to find sympy and/or numpy, ${PYTHON_SCRIPT} not running")
        return()
    endif()


    add_custom_command(
        OUTPUT
                ${PROJECT_SOURCE_DIR}/src/autogen/${OUTPUT_BASE}${OUTPUT_SUFFIX}
        COMMAND
                ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/autogen/${PYTHON_SCRIPT} ${PROJECT_SOURCE_DIR}/src/autogen/
        DEPENDS
                ${PROJECT_SOURCE_DIR}/src/autogen/${PYTHON_SCRIPT}
                ${PROJECT_SOURCE_DIR}/src/autogen/${TEMPLATE}
    )

    add_custom_target(
        autogen_${OUTPUT_BASE}
        DEPENDS
                ${PROJECT_SOURCE_DIR}/src/autogen/${OUTPUT_BASE}${OUTPUT_SUFFIX}
    )

    add_dependencies(${MAIN_TARGET} autogen_${OUTPUT_BASE})
endfunction()
