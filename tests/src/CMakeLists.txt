CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(OPESCI_CODEGEN)

# The version number.
set (OPESCI_VERSION_MAJOR 0)
set (OPESCI_VERSION_MINOR 1)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_language(CXX)

set (OPESCI_LIBRARIES)
# path of VTK in lab, ignored otherwise
set (VTK_DIR "/vol/bitbucket/ts2914/projects/VTK-build/lib/cmake/vtk-6.2")

FIND_PACKAGE(VTK NO_MODULE)

if(VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")
  include(${VTK_USE_FILE})
  set (OPESCI_LIBRARIES ${VTK_LIBRARIES} ${OPESCI_LIBRARIES})
  add_definitions(-DVTK_FOUND)
endif()
 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()

include_directories(../../include)

file(GLOB CXX_SOURCES ../../src/*.cpp)

# target_link_library
add_library(opesci SHARED ${CXX_SOURCES})
target_link_libraries(opesci ${OPESCI_LIBRARIES})

# add_subdirectory(tests EXCLUDE_FROM_ALL)

set(OPESCI_TEST_SRC ${CMAKE_SOURCE_DIR})

add_executable(test_grid3d ${OPESCI_TEST_SRC}/test_grid3d.cpp)
target_link_libraries(test_grid3d ${OPESCI_LIBRARIES} opesci)
add_executable(test_grid3d_read ${OPESCI_TEST_SRC}/test_grid3d_read.cpp)
target_link_libraries(test_grid3d_read ${OPESCI_LIBRARIES} opesci)
add_executable(test_grid3d_vtk ${OPESCI_TEST_SRC}/test_grid3d_vtk.cpp)
target_link_libraries(test_grid3d_vtk ${OPESCI_LIBRARIES} opesci)
