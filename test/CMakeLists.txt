itk_module_test()
add_executable(SkullStripping SkullStripping.cxx)
target_link_libraries(SkullStripping ${ITKSkullStrip_LIBRARIES})

set(TEMP ${ITK_TEST_OUTPUT_DIR})
itk_add_test( NAME SkullStrippingTest
  COMMAND SkullStripping
    DATA{Input/brainweb.mha}
    DATA{Input/atlasImage.mha}
    DATA{Input/atlasMask.mha}
    ${TEMP}/outputMask.mha
    ${TEMP}/outputMaskedImage.mha
)
