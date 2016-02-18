cmake_minimum_required(VERSION 2.8)
project(SM_standalone)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

include(ExternalProject)

if ( APPLE )
  set(CMAKE_C_FLAGS "-Wa,-q")
  set(CMAKE_CXX_FLAGS "-Wa,-q")
endif()

if ( WIN32 )
  add_compile_options(-g -O2)
else()
  add_compile_options(-g -fPIC -O2)
endif()

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
  BUILD_COMMAND make blosc_static
  CMAKE_ARGS -DBUILD_TESTS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
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

find_library(FFTW_LIBRARY NAMES fftw3 fftw fftw-3.3)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are strongly discouraged. "
    "Please remove CMakeCache.txt and the CMakeFiles/ directory, then build out-of-source.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -pthread")

if ( WIN32 )
  # necessary for FFTW to work on 32-bit, as mentioned e.g. here: https://github.com/libgit2/libgit2/issues/3342
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mincoming-stack-boundary=2") 
endif()

if ( WIN32 )
  # link all msys2/mingw libs statically on Windows
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic -static-libstdc++ -static-libgcc")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-Bstatic -static-libstdc++ -static-libgcc")

  # necessary for static compilation
  add_definitions(-DIN_LIBXML)
elseif ( APPLE )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
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
if (NOT APPLE)
  target_link_libraries(ims ${FFTW_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${BLOSC_LIBRARY})
endif()

# FIXME: get rid of boost
# on Mac OS X, homebrew hosts only clang-compiled binaries, but we'd rather use gcc
set(LIB_COMMON
  ims
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(convert tools/convert.cpp
  imzml/reader.cpp)
if ( WIN32 )
  target_link_libraries(convert -Wl,-Bstatic -lxml2 -lfftw3 -liconv ims -lws2_32)
elseif ( APPLE )
  target_link_libraries(convert -lxml2 -lfftw3 ims -lsnappy)
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
elseif( APPLE )
  target_link_libraries(detect ${Boost_PROGRAM_OPTIONS_LIBRARY} ims -lsnappy)
else()
  target_link_libraries(detect ${LIB_COMMON})
endif()
add_dependencies(detect blosc msgpack ims)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

add_library(ims_cffi SHARED cffi/image_measures.cpp cffi/imzb.cpp cffi/isocalc.cpp)
if ( WIN32 )
  target_link_libraries(ims_cffi -Wl,-Bstatic -lxml2 -lfftw3 -liconv ims -lwinpthread -lgomp -lws2_32)
elseif ( APPLE )
  target_link_libraries(ims_cffi ${BLOSC_LIBRARY} ims /usr/local/opt/fftw/lib/libfftw3.a)
else()
  target_link_libraries(ims_cffi ims)
endif()
add_dependencies(ims_cffi blosc ims)

# export compilation database for autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
