set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/bin)
set(OPESCI_TEST_SRC ${CMAKE_SOURCE_DIR}/tests/src)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests/include)

if(POCHOIR)
  set(LINK_LIBS "-L${CMAKE_BINARY_DIR}/lib;-lopesci")
  foreach(LIB ${OPESCI_LIBRARIES})
    set(LINK_LIBS "${LINK_LIBS};-l${LIB}")
  endforeach(LIB)

  set(pochoir_args -o bin/test_pochoir_iso_elastic -I${CMAKE_SOURCE_DIR}/include ${OPESCI_TEST_SRC}/test_pochoir_iso_elastic.cpp ${LINK_LIBS})
  add_custom_target(test_pochoir_iso_elastic ALL
                    COMMAND ${CMAKE_SOURCE_DIR}/tools/pochoir_compiler.sh ${pochoir_args}
                    SOURCES ${OPESCI_TEST_SRC}/test_pochoir_iso_elastic.cpp
                    VERBATIM)
  add_dependencies(test_pochoir_iso_elastic ${OPESCI_TEST_SRC}/test_pochoir_iso_elastic.cpp)

  add_test(test_pochoir_iso_elastic ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pochoir_iso_elastic)

  set_tests_properties(test_pochoir_iso_elastic PROPERTIES 
                       ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib")
endif()

add_executable(test_ref_iso_elastic ${OPESCI_TEST_SRC}/test_ref_iso_elastic.cpp)
target_link_libraries(test_ref_iso_elastic ${OPESCI_LIBRARIES} opesci)

add_test(test_ref_iso_elastic ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_ref_iso_elastic)
add_test(test_segy2vts ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_segy2vts)

