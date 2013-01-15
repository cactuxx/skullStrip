cmake_minimum_required(VERSION 2.8)


PROJECT(SkullStripping)


FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR"ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


add_executable(skullStripping skullStripping.cxx)
target_link_libraries(skullStripping ${ITK_LIBRARIES})