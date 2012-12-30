itk_module_test()
set(ITKSkullStripTests
SkullStripping.cxx
)

CreateTestDriver(ITKSkullStrip  "${ITKSkullStrip-Test_LIBRARIES}" "${ITKSkullStripTests}")

set(TEMP ${ITK_TEST_OUTPUT_DIR})
itk_add_test( NAME SkullStrippingTest
  COMMAND ITKSkullStripTestDriver
    --compare DATA{Baseline/outputMaskBaseline.mha}
              ${TEMP}/outputMask.mha
    SkullStripping
    DATA{Input/brainweb.mha}
    DATA{Input/atlasImage.mha}
    DATA{Input/atlasMask.mha}
    ${TEMP}/outputMask.mha
    ${TEMP}/outputMaskedImage.mha
)
