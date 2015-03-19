set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/bin)
set(OPESCI_TEST_SRC ${CMAKE_SOURCE_DIR}/tests/src)

#add_custom_target(test
#make
#COMMAND mkdir -p data
#COMMAND cp ${CMAKE_SOURCE_DIR}/tests/data/* data/
#COMMAND python ${CMAKE_SOURCE_DIR}/tests/unittest bin)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests/include)

ADD_EXECUTABLE(test_ref_iso_elastic ${OPESCI_TEST_SRC}/test_ref_iso_elastic.cpp)
TARGET_LINK_LIBRARIES(test_ref_iso_elastic opesci ${OPESCI_LIBRARIES})
