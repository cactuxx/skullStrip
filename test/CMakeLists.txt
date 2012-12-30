itk_module_test()
set(ITKSkullStripTests
SkullStripping.cxx
)

CreateTestDriver(ITKSkullStrip  "${ITKSkullStrip-Test_LIBRARIES}" "${ITKSkullStripTests}")

set(TEMP ${ITK_TEST_OUTPUT_DIR})
set(TEST_DATA_ROOT "${ITKSkullStrip_SOURCE_DIR}/test")
itk_add_test( NAME SkullStrippingTest
  COMMAND ITKSkullStripTestDriver
    --compare ${TEST_DATA_ROOT}/Baseline/outputMaskBaseline.mha
              ${TEMP}/outputMask.mha
    SkullStripping
    ${TEST_DATA_ROOT}/Input/brainweb.mha
    ${TEST_DATA_ROOT}/Input/atlasImage.mha
    ${TEST_DATA_ROOT}/Input/atlasMask.mha
    ${TEMP}/outputMask.mha
    ${TEMP}/outputMaskedImage.mha
)
