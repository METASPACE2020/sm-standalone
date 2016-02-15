cmake_minimum_required(VERSION 2.8)
project(SM_standalone)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

include(ExternalProject)

add_compile_options(-g -O2)

# MessagePack - simple serialization library
externalproject_add(msgpack
  GIT_REPOSITORY https://github.com/msgpack/msgpack-c
  GIT_TAG cpp-1.4.0
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  ) # for C++ compilation is not required
ExternalProject_Get_Property(msgpack source_dir)
include_directories(${source_dir}/include)

# blosc - fast compressor for numeric data
externalproject_add(blosc
  GIT_REPOSITORY https://github.com/Blosc/c-blosc
  GIT_TAG v1.7.1
  INSTALL_COMMAND ""
  CMAKE_ARGS -DBUILD_TESTS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
  )

ExternalProject_Get_Property(blosc source_dir)
include_directories(${source_dir}/blosc/)

find_package(OpenMP REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})

ExternalProject_Get_Property(blosc binary_dir)
set(BLOSC_LIBRARY ${binary_dir}/blosc/libblosc.a ${ZLIB_LIBRARY})

find_library(FFTW_LIBRARY NAMES fftw3 fftw)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are strongly discouraged. "
    "Please remove CMakeCache.txt and the CMakeFiles/ directory, then build out-of-source.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

if ( WIN32 )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic -static-libstdc++ -static-libgcc")
  add_definitions(-DIN_LIBXML)
endif()

set(SRC_COMMON
  imzb/reader.cpp
  imzb/writer.cpp
  ims/image_measures.cpp
  ms/isocalc.cpp
  ms/isotope_pattern.cpp
  ms/periodic_table.cpp
  )

add_library(ims STATIC ${SRC_COMMON})
target_link_libraries(ims ${FFTW_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${BLOSC_LIBRARY})

set(LIB_COMMON
  ims
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(convert tools/convert.cpp
  imzml/reader.cpp)
if ( WIN32 )
  target_link_libraries(convert -Wl,-Bstatic -lxml2 -lfftw3 -liconv ims -lws2_32)
else()
  target_link_libraries(convert ${LIBXML2_LIBRARIES} ${LIB_COMMON})
endif()
add_dependencies(convert blosc ims)

add_executable(isocalc tools/isocalc.cpp)
if ( WIN32 )
  target_link_libraries(isocalc -Wl,-Bstatic ${Boost_PROGRAM_OPTIONS_LIBRARY} -Wl,-Bstatic ims -lgomp -lfftw3 -lwinpthread -lws2_32)
else()
  target_link_libraries(isocalc ${LIB_COMMON})
endif()
add_dependencies(isocalc blosc msgpack ims)

add_executable(detect tools/detect.cpp)
if ( WIN32 )
  target_link_libraries(detect -Wl,-Bstatic ${Boost_PROGRAM_OPTIONS_LIBRARY} -Wl,-Bstatic ims -lgomp -lwinpthread -lws2_32)
else()
  target_link_libraries(detect ${LIB_COMMON})
endif()
add_dependencies(detect blosc msgpack ims)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

# export compilation database for Vim autocompletion with YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
