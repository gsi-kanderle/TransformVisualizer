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

#ifndef __qSlicerTransformVisualizerModuleWidget_h
#define __qSlicerTransformVisualizerModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"
#include "qSlicerTransformVisualizerModuleExport.h"

class qSlicerTransformVisualizerModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_TransformVisualizer
class Q_SLICER_QTMODULES_TRANSFORMVISUALIZER_EXPORT qSlicerTransformVisualizerModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT
  QVTK_OBJECT
  
public:
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerTransformVisualizerModuleWidget(QWidget *parent=0);
  virtual ~qSlicerTransformVisualizerModuleWidget();
  
  virtual void enter();

public slots:
  virtual void setMRMLScene(vtkMRMLScene*);
  void onSceneImportedEvent();
  void setTransformVisualizerParametersNode(vtkMRMLNode*);
  void update();
  void updateLabels();

  void visualize();

protected slots:
  void onLogicModified();
  void updateGlyphSourceOptions(int);
  
  void inputChanged(vtkMRMLNode*);
  void referenceVolumeChanged(vtkMRMLNode*);
  void outputModelChanged(vtkMRMLNode*);
  
  // Parameters
  // Glyph Parameters
  void setGlyphPointMax(double);
  void setGlyphScale(double);
  void setGlyphThreshold(double, double);
  void setGlyphSeed(int);
  void setSeed();
  void setGlyphSourceOption(int);
  // Arrow Parameters
  void setGlyphArrowScaleDirectional(bool);
  void setGlyphArrowScaleIsotropic(bool);  
  void setGlyphArrowTipLength(double);
  void setGlyphArrowTipRadius(double);
  void setGlyphArrowShaftRadius(double);
  void setGlyphArrowResolution(double);
  // Cone Parameters
  void setGlyphConeScaleDirectional(bool);
  void setGlyphConeScaleIsotropic(bool);  
  void setGlyphConeHeight(double);
  void setGlyphConeRadius(double);
  void setGlyphConeResolution(double);
  // Sphere Parameters
  void setGlyphSphereResolution(double);
    
  // Grid Parameters
  void setGridScale(double);
  void setGridSpacingMM(double);
  
  // Block Parameters
  void setBlockScale(double);
  void setBlockDisplacementCheck(int);
    
  // Contour Parameters
  void setContourValues(QString);
  void setContourDecimation(double);

  // Glyph Slice Parameters
  void setGlyphSliceNode(vtkMRMLNode*);
  void setGlyphSlicePointMax(double);
  void setGlyphSliceThreshold(double, double);
  void setGlyphSliceScale(double);
  void setGlyphSliceSeed(int);
  void setSeed2();
  
  // Grid Slice Parameters
  void setGridSliceNode(vtkMRMLNode*);
  void setGridSliceScale(double);
  void setGridSliceSpacingMM(double);
  
protected:
  QScopedPointer<qSlicerTransformVisualizerModuleWidgetPrivate> d_ptr;
  
  virtual void setup();
  void onEnter();

private:
  Q_DECLARE_PRIVATE(qSlicerTransformVisualizerModuleWidget);
  Q_DISABLE_COPY(qSlicerTransformVisualizerModuleWidget);
};

#endif
