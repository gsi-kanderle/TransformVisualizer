/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  This file was originally developed by Franklin King, PerkLab, Queen's University
  and was supported through the Applied Cancer Research Unit program of Cancer Care
  Ontario with funds provided by the Ontario Ministry of Health and Long-Term Care

==============================================================================*/

#ifndef __vtkSlicerTransformVisualizerLogic_h
#define __vtkSlicerTransformVisualizerLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>

// STD includes
#include <cstdlib>

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkMatrix4x4.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>

#include "vtkSlicerTransformVisualizerModuleLogicExport.h"

class vtkMRMLTransformVisualizerNode;
class vtkMRMLVectorVolumeNode;

/// \ingroup Slicer_QtModules_TransformVisualizer
class VTK_SLICER_TRANSFORMVISUALIZER_MODULE_LOGIC_EXPORT vtkSlicerTransformVisualizerLogic : 
  public vtkSlicerModuleLogic
{
public:
  static vtkSlicerTransformVisualizerLogic *New();
  vtkTypeMacro(vtkSlicerTransformVisualizerLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  enum visualizationModes
  {
    VIS_MODE_GLYPH_3D = 1,
    VIS_MODE_GRID_3D,
    VIS_MODE_BLOCK_3D,
    VIS_MODE_CONTOUR_3D,
    VIS_MODE_GLYPH_2D,
    VIS_MODE_GRID_2D
  };
 
  enum glyphSources
  {
    ARROW_3D = 0,
    CONE_3D,
    SPHERE_3D,
  }; 
  
  /*!
   * TODO: Add description of the function itself and its arguments
   */
  void CreateVisualization(int visualizationMode);
  void InitializeOutputModelNode(vtkMRMLModelNode* outputModelNode);
  
  void GlyphVisualization(bool inputIsTransform, vtkPolyData*, int);
  void GlyphSliceVisualization(bool inputIsTransform, vtkPolyData*); 
  void GlyphPreprocessInput(bool inputIsTransform, vtkUnstructuredGrid*, int seed, int pointMax, double min, double max);
  
  void GridVisualization(vtkImageData*, vtkPolyData* output);
  void GridSliceVisualization(vtkImageData*, vtkPolyData* output, vtkSmartPointer<vtkMatrix4x4>);
  
  void BlockVisualization(bool inputIsTransform, vtkPolyData* output);
  
  void ContourVisualization(bool inputIsTransform, vtkPolyData* output);
  
public:
  void SetAndObserveTransformVisualizerNode(vtkMRMLTransformVisualizerNode *node);
  vtkGetObjectMacro(TransformVisualizerNode, vtkMRMLTransformVisualizerNode);

protected:
  vtkSlicerTransformVisualizerLogic();
  ~vtkSlicerTransformVisualizerLogic();

  virtual void RegisterNodes();

private:
  vtkSlicerTransformVisualizerLogic(const vtkSlicerTransformVisualizerLogic&);// Not implemented
  void operator=(const vtkSlicerTransformVisualizerLogic&);// Not implemented

protected:
  /// Parameter set MRML node
  vtkMRMLTransformVisualizerNode* TransformVisualizerNode;
  
};

#endif
