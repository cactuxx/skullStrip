project(SkullStripping)

cmake_minimum_required( VERSION 2.6 )

#include some macros from another file...
include(${CMAKE_SOURCE_DIR}/IJMacros.txt)

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

#Build project executables
add_executable(skullStripping skullStripping.cxx itkStripTsImageFilter.h itkStripTsImageFilter.hxx)
target_link_libraries(skullStripping ${ITK_LIBRARIES})

if( ${ITK_MAJOR_VERSION} MATCHES "3" )
  add_executable(${CurrentExe} ImageCompare.cxx)
  target_link_libraries(${CurrentExe} ${Libraries})

  #Add testing
  add_test( NAME skullStrippingTest
    COMMAND skullStripping
      ${CMAKE_SOURCE_DIR}/brainweb.mha
      ${CMAKE_SOURCE_DIR}/atlasImage.mha
      ${CMAKE_SOURCE_DIR}/atlasMask.mha
    )

  add_test( NAME skullStrippingCompare
    COMMAND ImageCompare
      ${CMAKE_SOURCE_DIR}/outputMask.mha
      ${CMAKE_SOURCE_DIR}/outputMaskBaseline.mha
  )
else()

endif()
  add_test( NAME skullStrippingTest
    COMMAND skullStripping
      ${CMAKE_SOURCE_DIR}/brainweb.mha
      ${CMAKE_SOURCE_DIR}/atlasImage.mha
      ${CMAKE_SOURCE_DIR}/atlasMask.mha
    )

include(CTest)
