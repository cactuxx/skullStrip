#SkullStrip ITKv4 remote module
#contact: Stefan Bauer, stefan.bauer at istb.unibe.ch
#about: performs automatic skull-stripping, see Insight Journal article http://hdl.handle.net/10380/3353
 
itk_fetch_module(SkullStrip
  "A class to perform automatic skull-stripping for neuroimage analysis."
  GIT_REPOSITORY https://github.com/cactuxx/skullStrip.git
  GIT_TAG master
  )
