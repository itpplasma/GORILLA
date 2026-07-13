file(REMOVE
    "${ORBIT_DIR}/e_tot.dat"
    "${ORBIT_DIR}/full_orbit_plot_rphiz.dat"
    "${ORBIT_DIR}/p_phi.dat"
    "${ORBIT_DIR}/run.log")

execute_process(
    COMMAND "${GORILLA_MAIN}"
    WORKING_DIRECTORY "${ORBIT_DIR}"
    RESULT_VARIABLE run_status
    OUTPUT_VARIABLE run_stdout
    ERROR_VARIABLE run_stderr)
file(WRITE "${ORBIT_DIR}/run.log" "${run_stdout}${run_stderr}")
if(NOT run_status EQUAL 0)
    message(FATAL_ERROR "JOREK snapshot orbit failed; see ${ORBIT_DIR}/run.log")
endif()

execute_process(
    COMMAND "${ORBIT_VALIDATOR}"
    WORKING_DIRECTORY "${ORBIT_DIR}"
    RESULT_VARIABLE validation_status
    OUTPUT_VARIABLE validation_stdout
    ERROR_VARIABLE validation_stderr)
if(NOT validation_status EQUAL 0)
    message(FATAL_ERROR
        "JOREK snapshot orbit validation failed: ${validation_stdout}${validation_stderr}")
endif()
message(STATUS "${validation_stdout}")
