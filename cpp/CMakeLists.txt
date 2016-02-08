cmake_minimum_required (VERSION 2.8)
project (SM_standalone)

find_package(OpenMP)

find_package(LibXml2 REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})

find_library(FFTW_LIBRARY NAMES fftw3 fftw)
find_library(BLOSC_LIBRARY NAMES blosc)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are strongly discouraged. "
    "Please remove CMakeCache.txt and the CMakeFiles/ directory, then build out-of-source.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_compile_options(-O2)

set(SRC_COMMON
  ims/image_measures.cpp 
  ms/isocalc.cpp
  ms/isotope_pattern.cpp
  ms/periodic_table.cpp
  )

set(LIB_COMMON
  ${LIBXML2_LIBRARIES}
  ${FFTW_LIBRARY}
  ${BLOSC_LIBRARY}
  )

add_executable(convert tools/convert.cpp
  imzml/reader.cpp ${SRC_COMMON})
target_link_libraries(convert ${LIB_COMMON})

add_library(imzb SHARED cffi/imzb_cffi.cpp)
target_link_libraries(convert ${LIBXML2_LIBRARIES} ${BLOSC_LIBRARY})
