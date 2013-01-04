itk_module_test()
set(ITKSkullStripTests
SkullStrippingTest.cxx
)

CreateTestDriver(ITKSkullStrip  "${ITKSkullStrip-Test_LIBRARIES}" "${ITKSkullStripTests}")

# Build the sample application
set(ITK_NO_IO_FACTORY_REGISTER_MANAGER "")
include(${ITK_SOURCE_DIR}/CMake/UseITK.cmake)
add_executable(SkullStripping MACOSX_BUNDLE SkullStripping.cxx )
target_link_libraries(SkullStripping ${ITKSkullStrip_LIBRARIES})

set(TEMP ${ITK_TEST_OUTPUT_DIR})
set(TEST_DATA_ROOT "${ITKSkullStrip_SOURCE_DIR}/test")
itk_add_test( NAME SkullStrippingTest
  COMMAND ITKSkullStripTestDriver
    --compare ${TEST_DATA_ROOT}/Baseline/outputMaskBaseline.mha
              ${TEMP}/outputMask.mha
    SkullStrippingTest
    ${TEST_DATA_ROOT}/Input/brainweb.mha
    ${TEST_DATA_ROOT}/Input/atlasImage.mha
    ${TEST_DATA_ROOT}/Input/atlasMask.mha
    ${TEMP}/outputMask.mha
    ${TEMP}/outputMaskedImage.mha
)
