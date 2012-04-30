set(DOCUMENTATION "Here, write the documentation of the module" )

itk_module( ITKSkullStrip
  DEPENDS
    module_dep
  TEST_DEPENDS
    ITKTestKernel
    ITKIO
  DESCRIPTION
    "${DOCUMENTATION}"
)
