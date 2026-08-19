# Run Doxygen to completion and fail if it emits any diagnostic.

if(NOT DEFINED SOURCE_DIR)
    message(FATAL_ERROR "SOURCE_DIR is required.")
endif()

set(DOXYGEN_OUTPUT_DIR "${SOURCE_DIR}/docs/doxygen")
set(DOXYGEN_WARNING_LOG
    "${SOURCE_DIR}/docs/doxygen-warnings.log")

file(REMOVE_RECURSE "${DOXYGEN_OUTPUT_DIR}")
file(REMOVE "${DOXYGEN_WARNING_LOG}")

execute_process(
    COMMAND doxygen "${SOURCE_DIR}/Doxyfile"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE DOXYGEN_RESULT
    ERROR_VARIABLE DOXYGEN_STDERR
)

if(NOT "${DOXYGEN_STDERR}" STREQUAL "")
    message("${DOXYGEN_STDERR}")
endif()

set(DOXYGEN_LOG_CONTENT "")
if(EXISTS "${DOXYGEN_WARNING_LOG}")
    file(READ "${DOXYGEN_WARNING_LOG}" DOXYGEN_LOG_CONTENT)
endif()

if(NOT "${DOXYGEN_LOG_CONTENT}" STREQUAL "")
    message("${DOXYGEN_LOG_CONTENT}")
endif()

if(NOT "${DOXYGEN_RESULT}" STREQUAL "0")
    message(FATAL_ERROR
        "Doxygen exited with status ${DOXYGEN_RESULT}.")
endif()

if(NOT "${DOXYGEN_STDERR}" STREQUAL "" OR
   NOT "${DOXYGEN_LOG_CONTENT}" STREQUAL "")
    message(FATAL_ERROR
        "Doxygen emitted warnings or errors.")
endif()
