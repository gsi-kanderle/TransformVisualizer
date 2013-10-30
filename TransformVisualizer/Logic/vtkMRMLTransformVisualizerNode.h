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

#pragma once
#ifndef __vtkMRMLTransformVisualizerNode_h
#define __vtkMRMLTransformVisualizerNode_h

#include <vtkMRML.h>
#include <vtkMRMLNode.h>

#include "vtkSlicerTransformVisualizerModuleLogicExport.h"

/// \ingroup Slicer_QtModules_TransformVisualizer
class VTK_SLICER_TRANSFORMVISUALIZER_MODULE_LOGIC_EXPORT vtkMRMLTransformVisualizerNode : 
  public vtkMRMLNode
{
public:   

  static vtkMRMLTransformVisualizerNode *New();
  vtkTypeMacro(vtkMRMLTransformVisualizerNode, vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  static std::string InputReferenceRole;
  static std::string ReferenceVolumeReferenceRole;
  static std::string OutputModelReferenceRole;

  static std::string GlyphSliceReferenceRole;
  static std::string GridSliceReferenceRole;

  virtual vtkMRMLNode* CreateNodeInstance();
  virtual void ReadXMLAttributes( const char** atts);
  virtual void WriteXML(ostream& of, int indent);
  virtual void Copy(vtkMRMLNode *node);
  virtual const char* GetNodeTagName() {return "TransformVisualizerParameters";};

public:
  vtkMRMLNode* GetInputNode();
  void SetAndObserveInputNode(vtkMRMLNode* node);

  vtkMRMLNode* GetReferenceVolumeNode();
  void SetAndObserveReferenceVolumeNode(vtkMRMLNode* node);

  vtkMRMLNode* GetOutputModelNode();
  void SetAndObserveOutputModelNode(vtkMRMLNode* node);

  // Glyph Parameters
  vtkSetMacro(GlyphPointMax, int);
  vtkGetMacro(GlyphPointMax, int);
  vtkSetMacro(GlyphScale, float);
  vtkGetMacro(GlyphScale, float);
  vtkSetMacro(GlyphThresholdMax, double);
  vtkGetMacro(GlyphThresholdMax, double);
  vtkSetMacro(GlyphThresholdMin, double);
  vtkGetMacro(GlyphThresholdMin, double);  
  vtkSetMacro(GlyphSeed, int);
  vtkGetMacro(GlyphSeed, int);
  vtkSetMacro(GlyphSourceOption, int);
  vtkGetMacro(GlyphSourceOption, int);
  // Arrow Parameters
  vtkSetMacro(GlyphArrowScaleDirectional, bool);
  vtkGetMacro(GlyphArrowScaleDirectional, bool);
  vtkSetMacro(GlyphArrowScaleIsotropic, bool);
  vtkGetMacro(GlyphArrowScaleIsotropic, bool);  
  vtkSetMacro(GlyphArrowTipLength, float);
  vtkGetMacro(GlyphArrowTipLength, float);
  vtkSetMacro(GlyphArrowTipRadius, float);
  vtkGetMacro(GlyphArrowTipRadius, float);
  vtkSetMacro(GlyphArrowShaftRadius, float);
  vtkGetMacro(GlyphArrowShaftRadius, float);
  vtkSetMacro(GlyphArrowResolution, int);
  vtkGetMacro(GlyphArrowResolution, int);
  // Cone Parameters
  vtkSetMacro(GlyphConeScaleDirectional, bool);
  vtkGetMacro(GlyphConeScaleDirectional, bool);
  vtkSetMacro(GlyphConeScaleIsotropic, bool);
  vtkGetMacro(GlyphConeScaleIsotropic, bool);    
  vtkSetMacro(GlyphConeHeight, float);
  vtkGetMacro(GlyphConeHeight, float);
  vtkSetMacro(GlyphConeRadius, float);
  vtkGetMacro(GlyphConeRadius, float);
  vtkSetMacro(GlyphConeResolution, int);
  vtkGetMacro(GlyphConeResolution, int);
  // Sphere Parameters
  vtkSetMacro(GlyphSphereResolution, float);
  vtkGetMacro(GlyphSphereResolution, float);
    
  // Grid Parameters
  vtkSetMacro(GridScale, float);
  vtkGetMacro(GridScale, float);
  vtkSetMacro(GridSpacingMM, int);
  vtkGetMacro(GridSpacingMM, int);
  
  // Block Parameters
  vtkSetMacro(BlockScale, float);
  vtkGetMacro(BlockScale, float);
  vtkSetMacro(BlockDisplacementCheck, int);
  vtkGetMacro(BlockDisplacementCheck, int);  
    
  // Contour Parameters
  vtkSetMacro(ContourNumber, int);
  vtkGetMacro(ContourNumber, int);
  void SetContourValues(double*, int size);
  double* GetContourValues();
  vtkSetMacro(ContourDecimation, float);
  vtkGetMacro(ContourDecimation, float); 
  
  // Glyph Slice Parameters
  vtkMRMLNode* GetGlyphSliceNode();
  void SetAndObserveGlyphSliceNode(vtkMRMLNode* node);
  vtkSetMacro(GlyphSlicePointMax, int);
  vtkGetMacro(GlyphSlicePointMax, int);
  vtkSetMacro(GlyphSliceThresholdMax, double);
  vtkGetMacro(GlyphSliceThresholdMax, double);
  vtkSetMacro(GlyphSliceThresholdMin, double);
  vtkGetMacro(GlyphSliceThresholdMin, double);    
  vtkSetMacro(GlyphSliceScale, float);
  vtkGetMacro(GlyphSliceScale, float);
  vtkSetMacro(GlyphSliceSeed, int);
  vtkGetMacro(GlyphSliceSeed, int);
  
  //Grid Slice Parameters
  vtkMRMLNode* GetGridSliceNode();
  void SetAndObserveGridSliceNode(vtkMRMLNode* node);
  vtkSetMacro(GridSliceScale, float);
  vtkGetMacro(GridSliceScale, float);
  vtkSetMacro(GridSliceSpacingMM, int);
  vtkGetMacro(GridSliceSpacingMM, int);
  
protected:
  vtkMRMLTransformVisualizerNode();
  ~vtkMRMLTransformVisualizerNode();

  vtkMRMLTransformVisualizerNode(const vtkMRMLTransformVisualizerNode&);
  void operator=(const vtkMRMLTransformVisualizerNode&);

//Parameters
protected:
  // Glyph Parameters
  int GlyphPointMax;
  //TODO: Need to change the UI into float too
  float GlyphScale;
  double GlyphThresholdMax;
  double GlyphThresholdMin;
  int GlyphSeed;
  int GlyphSourceOption; //0 - Arrow, 1 - Cone, 2 - Sphere
  // Arrow Parameters
  bool GlyphArrowScaleDirectional;
  bool GlyphArrowScaleIsotropic;  
  float GlyphArrowTipLength;
  float GlyphArrowTipRadius;
  float GlyphArrowShaftRadius;
  int GlyphArrowResolution;
  
  // Cone Parameters
  bool GlyphConeScaleDirectional;
  bool GlyphConeScaleIsotropic;    
  float GlyphConeHeight;
  float GlyphConeRadius;
  int GlyphConeResolution;
  
  // Sphere Parameters
  float GlyphSphereResolution;

  // Grid Parameters
  float GridScale;
  int GridSpacingMM;
  
  // Block Parameters
  float BlockScale;
  int BlockDisplacementCheck;
    
  // Contour Parameters
  int ContourNumber;
  double* ContourValues;
  float ContourDecimation;

  // Glyph Slice Parameters
  int GlyphSlicePointMax;
  double GlyphSliceThresholdMax;
  double GlyphSliceThresholdMin;
  float GlyphSliceScale;
  int GlyphSliceSeed;
  
  // Grid Slice Parameters
  float GridSliceScale;
  int GridSliceSpacingMM;
};

#endif
