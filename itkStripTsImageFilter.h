/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#ifndef __itkStripTsImageFilter_h
#define __itkStripTsImageFilter_h


#include "itkImageToImageFilter.h"

#include "itkImageDuplicator.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkCastImageFilter.h"


namespace itk
{

/** \class StripTsImageFilter
 * \brief composite ITK filter for automatic skull-stripping
 *
 * This ITK filter performs automatic skull stripping for 3D medical images
 *
 * It requires 3 inputs:
 * SetInput()
 * SetAtlasImage()
 * SetAtlasBrainMask()
 * and it outputs the brain mask for the patient image
 *
 * \warning images have to be 3D
 *
 *
 * to be done: check if all required inputs are set and if images are 3D
 *
 * version: 1.0
 * date: 20.04.2012
 * initial release

 * Stefan Bauer
 * Institute for Surgical Technology and Biomechanics, University of Bern
 * stefan.bauer@istb.unibe.ch
 *
 */

template <class TImageType, class TAtlasImageType, class TAtlasLabelType>
class ITK_EXPORT
StripTsImageFilter : public ImageToImageFilter<TImageType, TAtlasLabelType>
{
public:

  // standard class typedefs
  typedef StripTsImageFilter                              Self;
  typedef ImageToImageFilter<TImageType,TAtlasLabelType>  Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  // method for creation through the object factory
  itkNewMacro(Self);

  // run-time type information (and related methods)
  itkTypeMacro(StripTsImageFilter, ImageToImageFilter);

  // display
  void PrintSelf( std::ostream& os, Indent indent ) const;

  // image and label templates
  typedef TImageType                        ImageType;
  typedef typename ImageType::Pointer       ImagePointer;
  typedef typename ImageType::ConstPointer  ImageConstPointer;

  typedef TAtlasImageType                       AtlasImageType;
  typedef typename AtlasImageType::Pointer      AtlasImagePointer;
  typedef typename AtlasImageType::ConstPointer AtlasImageConstPointer;

  typedef TAtlasLabelType                       AtlasLabelType;
  typedef typename AtlasLabelType::Pointer      AtlasLabelPointer;
  typedef typename AtlasLabelType::ConstPointer AtlasLabelConstPointer;

  void SetAtlasImage( const TAtlasImageType * ptr );
  void SetAtlasBrainMask( const TAtlasLabelType * ptr );


protected:

  StripTsImageFilter();
  ~StripTsImageFilter(){}

  // does the real work
  virtual void GenerateData();


private:

  StripTsImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ImagePointer      m_PatientImage;
  AtlasImagePointer m_AtlasImage;
  AtlasLabelPointer m_AtlasLabels;

  void RescaleImages();
  void DownsampleImage();
  void RigidRegistration();
  void AffineRegistration();
  void BinaryErosion();
  void MultiResLevelSet();
  void PyramidFilter(int isoSpacing);
  void InversePyramidFilter();
  void LevelSetRefinement(int isoSpacing);
  void UpsampleLabels();

}; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStripTsImageFilter.hxx"
#endif

#endif
