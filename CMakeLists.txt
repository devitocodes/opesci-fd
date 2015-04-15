CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(OPESCI)

# The version number.
set (OPESCI_VERSION_MAJOR 0)
set (OPESCI_VERSION_MINOR 1)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(CTest)
enable_testing()

enable_language(CXX)

find_program(POCHOIR pochoir PATHS)
if(POCHOIR)
  MESSAGE(STATUS "Found the Pochoir Stencil Compiler: ${POCHOIR}")
else()
  MESSAGE(WARNING "WARNING: Cannot find the Pochoir Stencil Compiler. Pochoir source will not be built.")
endif()

set (OPESCI_LIBRARIES)
FIND_PACKAGE(VTK NO_MODULE)
if(VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")

  include(${VTK_USE_FILE})
  set (OPESCI_LIBRARIES ${VTK_LIBRARIES} ${OPESCI_LIBRARIES})
  add_definitions(-DVTK_FOUND)
endif()
 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-deprecated")
endif()
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  message(STATUS "Adding magic sauce for Intel compiler.")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -xHOST -g")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()

include_directories(include)

file(GLOB CXX_SOURCES src/*.cpp)

# target_link_library
add_library(opesci SHARED ${CXX_SOURCES})
target_link_libraries(opesci)

add_subdirectory(tests EXCLUDE_FROM_ALL)

install(DIRECTORY include/ DESTINATION include/opesci FILES_MATCHING PATTERN *.h)
install(FILES libopesci.so DESTINATION lib)


