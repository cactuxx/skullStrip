PROJECT(SkullStripping)

#include some macros from another file...
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)


#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)


#Declare any external dependencies that your project may have here.
SET(Required_Packages
ITK
)


#this foreach loads all of the packages that you specified as required.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)


#Set any libraries that your project depends on.
SET(Libraries
${ITK_LIBRARIES}
)


#Build project executables
SET(CurrentExe "skullStripping")
ADD_EXECUTABLE(${CurrentExe} skullStripping.cxx itkStripTsImageFilter.h itkStripTsImageFilter.hxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "imageCompare")
ADD_EXECUTABLE(${CurrentExe} imageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})


#Add testing
ADD_TEST(skullStrippingTest skullStripping ${CMAKE_SOURCE_DIR}/brainweb.mha ${CMAKE_SOURCE_DIR}/atlasImage.mha ${CMAKE_SOURCE_DIR}/atlasMask.mha)
ADD_TEST(skullStrippingCompare imageCompare ${CMAKE_SOURCE_DIR}/outputMask.mha ${CMAKE_SOURCE_DIR}/outputMaskBaseline.mha)
