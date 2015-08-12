set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/bin)
set(OPESCI_TEST_SRC ${CMAKE_SOURCE_DIR}/tests/src)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests/include)

add_executable(test_ref_iso_elastic ${OPESCI_TEST_SRC}/test_ref_iso_elastic.cpp)
target_link_libraries(test_ref_iso_elastic ${OPESCI_LIBRARIES} opesci)

add_executable(test1 ${OPESCI_TEST_SRC}/test1.cpp)
target_link_libraries(test1 ${OPESCI_LIBRARIES} opesci)


add_executable(test_segy2vts ${OPESCI_TEST_SRC}/test_segy2vts.cpp)
target_link_libraries(test_segy2vts ${OPESCI_LIBRARIES} opesci)

add_test(test_ref_iso_elastic ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_ref_iso_elastic)
add_test(test_segy2vts ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_segy2vts)

