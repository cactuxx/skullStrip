project(SkullStripping)

cmake_minimum_required( VERSION 2.6 )

#include some macros from another file...
include(${CMAKE_SOURCE_DIR}/IJMacros.txt)

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

include(CTest)

include_directories( ${CMAKE_SOURCE_DIR}/include )

add_subdirectory( test )
