set(DOCUMENTATION "This module contains a class to perform automatic skull-stripping for neuroimage analysis.
It is based on the ITK level-set and registration frameworks." )

itk_module( ITKSkullStrip
  DEPENDS
    ITKLevelSets 
    ITKRegistrationCommon
  TEST_DEPENDS
    ITKTestKernel
    ITKIOMeta
  DESCRIPTION
    "${DOCUMENTATION}"
)
