function(prepend_path PATH SOURCE_FILES)
        # Use recursive substitution to expand SOURCE_FILES
        unset(MODIFIED)
        foreach(SOURCE_FILE IN ITEMS ${${SOURCE_FILES}})
                list(APPEND MODIFIED "${PATH}/${SOURCE_FILE}")
        endforeach()
        set(${SOURCE_FILES} ${MODIFIED} PARENT_SCOPE)
endfunction()
