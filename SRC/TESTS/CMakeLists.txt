set(CMAKE_VERBOSE_MAKEFILE ON)

add_pfunit_ctest (test_various_functions_mod
  TEST_SOURCES test_various_functions_mod.pf
  LINK_LIBRARIES GORILLA
  )

add_pfunit_ctest (test_binsrc
  TEST_SOURCES test_binsrc.pf
  LINK_LIBRARIES GORILLA
  )
