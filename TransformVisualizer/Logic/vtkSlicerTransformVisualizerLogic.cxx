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

// TransformVisualizer includes
#include "vtkMRMLTransformVisualizerNode.h"
#include "vtkSlicerTransformVisualizerLogic.h"
#include "vtkTransformVisualizerGlyph3D.h"

// MRML includes
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLColorTableNode.h>
#include <vtkMRMLTransformNode.h>
#include <vtkMRMLSliceNode.h>

// VTK includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkVectorNorm.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkGeneralTransform.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>

// Glyph VTK includes
#include <vtkArrowSource.h>
#include <vtkConeSource.h>
#include <vtkSphereSource.h>

// Grid VTK includes
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkWarpVector.h>
#include <vtkExtractVOI.h>
#include <vtkImageResample.h>

// Block VTK includes
#include <vtkGeometryFilter.h>
#include <vtkTransformFilter.h>
#include <vtkVectorDot.h>
#include <vtkPolyDataNormals.h>

// Contour VTK includes
#include <vtkMarchingCubes.h>
#include <vtkDecimatePro.h>

// Glyph Slice VTK includes
#include <vtkGlyphSource2D.h>
#include <vtkRibbonFilter.h>
#include <vtkPlane.h>

// Grid Slice VTK includes
#include <vtkAppendPolyData.h>

// STD includes
#include <cassert>
#include <math.h>

#include <vtkTimerLog.h>

#define EPSILON 0.0001

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerTransformVisualizerLogic);

//----------------------------------------------------------------------------
vtkSlicerTransformVisualizerLogic::vtkSlicerTransformVisualizerLogic()
{
  this->TransformVisualizerNode = NULL;
}

//----------------------------------------------------------------------------
vtkSlicerTransformVisualizerLogic::~vtkSlicerTransformVisualizerLogic()
{
  vtkSetAndObserveMRMLNodeMacro(this->TransformVisualizerNode, NULL);
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::SetAndObserveTransformVisualizerNode(vtkMRMLTransformVisualizerNode *node)
{
  vtkSetAndObserveMRMLNodeMacro(this->TransformVisualizerNode, node);
}

//-----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::RegisterNodes()
{
  vtkMRMLScene* scene = this->GetMRMLScene();
  if (scene == NULL)
  {
    vtkErrorMacro("Null scene");
    return;
  }
  
  scene->RegisterNodeClass(vtkSmartPointer<vtkMRMLTransformVisualizerNode>::New());
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::CreateVisualization(int visualizationMode)
{
  if (this->TransformVisualizerNode == NULL || this->GetMRMLScene() == NULL)
  {
    vtkErrorMacro("CreateVisualization failed: Deformation Field Visualizer Node or scene is invalid");
    return;
  }
 
  vtkMRMLModelNode* outputModelNode = vtkMRMLModelNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetOutputModelNodeID()));
  if (outputModelNode == NULL)
  {
    vtkErrorMacro("CreateVisualization failed: Model node is invalid");
    return;
  }
  
  bool inputIsTransform;
  vtkMRMLNode* initialInputNode = this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID());
  
  if (initialInputNode == NULL)
  {
    vtkErrorMacro("CreateVisualization failed: Input node is invalid");
    return;
  }
  if (initialInputNode->IsA("vtkMRMLLinearTransformNode") || initialInputNode->IsA("vtkMRMLBSplineTransformNode") || initialInputNode->IsA("vtkMRMLGridTransformNode"))
  {
    if (this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetReferenceVolumeNodeID()) == NULL)
    {
      vtkErrorMacro("CreateVisualization failed: Reference volume node is invalid");
      return;
    }
    inputIsTransform = true;
  }  
  else if (initialInputNode->IsA("vtkMRMLVectorVolumeNode"))
  {
    inputIsTransform = false;
  }
  else
  {
    vtkErrorMacro("Invalid input node selected. Expected vtkMRMLVectorVolumeNode, vtkMRMLLinearTransformNode, vtkMRMLBSplineTransformNode, or vtkMRMLGridTransformNode, but got" << (this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()))->GetClassName() << "instead");
    return;
  }
  
  this->GetMRMLScene()->StartState(vtkMRMLScene::BatchProcessState); 
  
  this->InitializeOutputModelNode(outputModelNode);

  vtkPolyData* output = vtkPolyData::New();
  switch (visualizationMode){
    case VIS_MODE_GLYPH_3D:
      this->GlyphVisualization(inputIsTransform, output, this->TransformVisualizerNode->GetGlyphSourceOption());
      outputModelNode->GetModelDisplayNode()->SetScalarVisibility(1);
      outputModelNode->GetModelDisplayNode()->SetActiveScalarName("VectorMagnitude");
      break;
    case VIS_MODE_GLYPH_2D:
      if (vtkMRMLSliceNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetGlyphSliceNodeID())) == NULL)
      {
        vtkErrorMacro("Failed to create Glyph Slice visualization: Invalid slice node");
        return;
      }
      this->GlyphSliceVisualization(inputIsTransform, output);
      outputModelNode->GetModelDisplayNode()->SetScalarVisibility(1);
      outputModelNode->GetModelDisplayNode()->SetActiveScalarName("OriginalVectorMagnitude");
      outputModelNode->GetModelDisplayNode()->SetBackfaceCulling(0);
      outputModelNode->GetModelDisplayNode()->SetSliceIntersectionVisibility(1);
      break;     
    case VIS_MODE_BLOCK_3D:
      this->BlockVisualization(inputIsTransform, output);
      outputModelNode->GetModelDisplayNode()->SetScalarVisibility(1);
      outputModelNode->GetModelDisplayNode()->SetActiveScalarName("VectorMagnitude");
      outputModelNode->GetModelDisplayNode()->SetBackfaceCulling(0);
      break;
  }
  
  outputModelNode->SetAndObservePolyData(output);
  output->Delete();
  
  this->GetMRMLScene()->EndState(vtkMRMLScene::BatchProcessState);
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::InitializeOutputModelNode(vtkMRMLModelNode* outputModelNode)
{
  // Output Node
  if (outputModelNode->GetModelDisplayNode()==NULL)
  {
    vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode = vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    displayNode = vtkMRMLModelDisplayNode::SafeDownCast(this->GetMRMLScene()->AddNode(displayNode));
    outputModelNode->SetAndObserveDisplayNodeID(displayNode->GetID());
  }  
  
  outputModelNode->SetHideFromEditors(0);
  outputModelNode->SetSelectable(1);
  outputModelNode->Modified();
  
  if (outputModelNode->GetModelDisplayNode()->GetColorNode()==NULL)
  {
    vtkSmartPointer<vtkMRMLColorTableNode> colorTableNode = vtkSmartPointer<vtkMRMLColorTableNode>::New();
    this->GetMRMLScene()->AddNode(colorTableNode);
    
    colorTableNode->SetName("Deformation Field Colors");
    colorTableNode->SetAttribute("Category", "User Generated");
    colorTableNode->SetTypeToUser();
    colorTableNode->SetNumberOfColors(4);
    colorTableNode->GetLookupTable();
    colorTableNode->AddColor("negligible", 0.0, 0.0, 0.5, 1.0);
    colorTableNode->AddColor(       "low", 0.0, 1.0, 0.0, 1.0);
    colorTableNode->AddColor(    "medium", 1.0, 1.0, 0.0, 1.0);
    colorTableNode->AddColor(      "high", 1.0, 0.0, 0.0, 1.0);

    outputModelNode->GetModelDisplayNode()->SetAndObserveColorNodeID(colorTableNode->GetID());
  }
  
  vtkMRMLColorTableNode *colorNode = vtkMRMLColorTableNode::SafeDownCast(outputModelNode->GetModelDisplayNode()->GetColorNode());
  /*
  if (colorNode != NULL)
  {
	//How to do after new input process?
  //For now just do with what there is
    double* range = deformationField->GetPointData()->GetScalars()->GetRange(-1);
    colorNode->GetLookupTable()->SetTableRange(range[0], range[1]);
  }
  */
}

//Glyph Visualization
//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::GlyphVisualization(bool inputIsTransform, vtkPolyData* output, int sourceOption)
{
  //Pre-processing
  vtkSmartPointer<vtkUnstructuredGrid> pointSet = vtkSmartPointer<vtkUnstructuredGrid>::New();
  pointSet->Initialize();
  this->GlyphPreprocessInput(inputIsTransform, pointSet, this->TransformVisualizerNode->GetGlyphSeed(), this->TransformVisualizerNode->GetGlyphPointMax(),
                            this->TransformVisualizerNode->GetGlyphThresholdMin(), this->TransformVisualizerNode->GetGlyphThresholdMax());
  
  vtkSmartPointer<vtkTransformVisualizerGlyph3D> glyphFilter = vtkSmartPointer<vtkTransformVisualizerGlyph3D>::New();
  glyphFilter->SetScaleModeToScaleByVector();
  glyphFilter->SetScaleFactor(this->TransformVisualizerNode->GetGlyphScale());
  glyphFilter->SetColorModeToColorByVector();

  switch (sourceOption){
    //Arrows
    case ARROW_3D:
    {
      glyphFilter->SetScaleDirectional(this->TransformVisualizerNode->GetGlyphArrowScaleDirectional());
      vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
      arrowSource->SetTipLength(this->TransformVisualizerNode->GetGlyphArrowTipLength());
      arrowSource->SetTipRadius(this->TransformVisualizerNode->GetGlyphArrowTipRadius());
      arrowSource->SetTipResolution(this->TransformVisualizerNode->GetGlyphArrowResolution());
      arrowSource->SetShaftRadius(this->TransformVisualizerNode->GetGlyphArrowShaftRadius());
      arrowSource->SetShaftResolution(this->TransformVisualizerNode->GetGlyphArrowResolution());
      
      glyphFilter->OrientOn();
      glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
      break;
    }
    //Cones
    case CONE_3D:
    {
      glyphFilter->SetScaleDirectional(this->TransformVisualizerNode->GetGlyphConeScaleDirectional());
      vtkSmartPointer<vtkConeSource> coneSource = vtkSmartPointer<vtkConeSource>::New();
      coneSource->SetHeight(this->TransformVisualizerNode->GetGlyphConeHeight());
      coneSource->SetRadius(this->TransformVisualizerNode->GetGlyphConeRadius());
      coneSource->SetResolution(this->TransformVisualizerNode->GetGlyphConeResolution());
      
      glyphFilter->OrientOn();
      glyphFilter->SetSourceConnection(coneSource->GetOutputPort());
      break;
    }
    //Spheres
    case SPHERE_3D:
    {
      glyphFilter->SetScaleDirectional(false);
      vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource->SetRadius(1);
      sphereSource->SetThetaResolution(this->TransformVisualizerNode->GetGlyphSphereResolution());
      sphereSource->SetPhiResolution(this->TransformVisualizerNode->GetGlyphSphereResolution());
      
      glyphFilter->OrientOn();
      glyphFilter->SetSourceConnection(sphereSource->GetOutputPort());
      break;
    }
  }
  glyphFilter->SetInputConnection(pointSet->GetProducerPort());
  glyphFilter->Update();
  
  output->ShallowCopy(glyphFilter->GetOutput());
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::GlyphSliceVisualization(bool inputIsTransform, vtkPolyData* output)
{
  //Pre-processing
  vtkSmartPointer<vtkUnstructuredGrid> pointSet = vtkSmartPointer<vtkUnstructuredGrid>::New();
  pointSet->Initialize();
  this->GlyphPreprocessInput(inputIsTransform, pointSet, this->TransformVisualizerNode->GetGlyphSliceSeed(), this->TransformVisualizerNode->GetGlyphSlicePointMax(), 
                            this->TransformVisualizerNode->GetGlyphSliceThresholdMin(), this->TransformVisualizerNode->GetGlyphSliceThresholdMax());
  

  vtkSmartPointer<vtkRibbonFilter> ribbon = vtkSmartPointer<vtkRibbonFilter>::New();
  float sliceNormal[3] = {0,0,0};
  double width = 1;
  const double* spacing = NULL;
  
  vtkMRMLSliceNode* sliceNode = vtkMRMLSliceNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetGlyphSliceNodeID()));
  vtkSmartPointer<vtkMatrix4x4> ijkToRasDirections = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkMatrix4x4> rasToIjkDirections = vtkSmartPointer<vtkMatrix4x4>::New();
  if (inputIsTransform)
  {
    vtkMRMLVolumeNode* referenceVolumeNode = vtkMRMLVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetReferenceVolumeNodeID()));
    referenceVolumeNode->GetIJKToRASDirectionMatrix(ijkToRasDirections);
    spacing = referenceVolumeNode->GetImageData()->GetSpacing();
  }
  else
  {
    vtkMRMLVectorVolumeNode* inputVectorVolumeNode = vtkMRMLVectorVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    inputVectorVolumeNode->GetIJKToRASDirectionMatrix(ijkToRasDirections);
    spacing = inputVectorVolumeNode->GetSpacing();
  }
  rasToIjkDirections->DeepCopy(ijkToRasDirections);
  vtkSmartPointer<vtkMatrix4x4> sliceToIjk = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMatrix4x4::Multiply4x4(rasToIjkDirections,sliceNode->GetSliceToRAS(),sliceToIjk);
  
  sliceNormal[0] = sliceToIjk->GetElement(0,2);
  sliceNormal[1] = sliceToIjk->GetElement(1,2);
  sliceNormal[2] = sliceToIjk->GetElement(2,2);
  
  if (abs(sliceNormal[0]) > EPSILON)
  {
    width = spacing[0];
  }
  else if (abs(sliceNormal[1]) > EPSILON)
  {
    width = spacing[1];
  }    
  else if (abs(sliceNormal[2]) > EPSILON)
  {
    width = spacing[2];
  }
  //Oblique
  //TODO: Currently takes the largest spacing to avoid any issues. Only affects oblique slices, but will change later to a more precise method when slice options are updated
  else
  {
    if (spacing[0] >= spacing[1] && spacing[0] >= spacing[2])
    {
      width = spacing[0];
    }
    else if (spacing[1] >= spacing[0] && spacing[1] >= spacing[2])
    {
      width = spacing[1];
    }
    else
    {
      width = spacing[2];
    }
  }

  vtkSmartPointer<vtkVectorNorm> norm = vtkSmartPointer<vtkVectorNorm>::New();
  norm->SetInputConnection(pointSet->GetProducerPort());
  norm->Update();
  vtkSmartPointer<vtkDataArray> vectorMagnitude = norm->GetOutput()->GetPointData()->GetScalars();  
  
  //Projection to slice plane
  float dot = 0;
  vtkDataArray* projectedVectors = pointSet->GetPointData()->GetVectors();
  double* chosenVector = NULL;
  for(int i = 0; i < projectedVectors->GetNumberOfTuples(); i++)
  {
    chosenVector = projectedVectors->GetTuple3(i);
    dot = chosenVector[0]*sliceNormal[0] + chosenVector[1]*sliceNormal[1] + chosenVector[2]*sliceNormal[2];
    projectedVectors->SetTuple3(i, chosenVector[0]-dot*sliceNormal[0], chosenVector[1]-dot*sliceNormal[1], chosenVector[2]-dot*sliceNormal[2]);
  }
  
  vectorMagnitude->SetName("OriginalVectorMagnitude");
  pointSet->GetPointData()->AddArray(vectorMagnitude);
    
  vtkSmartPointer<vtkTransform> rotateArrow = vtkSmartPointer<vtkTransform>::New();
  rotateArrow->RotateX(vtkMath::DegreesFromRadians(acos(abs(sliceNormal[2]))));
  
  vtkSmartPointer<vtkGlyphSource2D> arrow2DSource = vtkSmartPointer<vtkGlyphSource2D>::New();
  arrow2DSource->SetGlyphTypeToArrow();
  arrow2DSource->SetScale(1);
  arrow2DSource->SetFilled(0);
  
  vtkSmartPointer<vtkTransformVisualizerGlyph3D> glyphFilter = vtkSmartPointer<vtkTransformVisualizerGlyph3D>::New();
  glyphFilter->SetScaleModeToScaleByVector();
  glyphFilter->SetScaleFactor(this->TransformVisualizerNode->GetGlyphSliceScale());
  glyphFilter->SetScaleDirectional(false);
  glyphFilter->SetColorModeToColorByVector();
  glyphFilter->SetSourceTransform(rotateArrow);
  
  glyphFilter->SetSourceConnection(arrow2DSource->GetOutputPort());
  glyphFilter->SetInputConnection(pointSet->GetProducerPort());
  glyphFilter->Update();
  
  ribbon->SetInputConnection(glyphFilter->GetOutputPort());
  ribbon->SetDefaultNormal(sliceNormal[0], sliceNormal[1], sliceNormal[2]);
  ribbon->SetWidth(width/2 + 0.15);
  ribbon->SetAngle(90.0);
  ribbon->UseDefaultNormalOn();
  ribbon->Update();
  
  output->ShallowCopy(ribbon->GetOutput());
}

//Preprocess Glyph
//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::GlyphPreprocessInput(bool inputIsTransform, vtkUnstructuredGrid* pointSet, int seed, int pointMax, double min, double max)
{
  //Will contain all the points that are to be rendered
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  
  //Will contain the corresponding vectors for pointSet
  vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
  vectors->Initialize();
  vectors->SetNumberOfComponents(3);  

  vtkSmartPointer<vtkMatrix4x4> IJKToRAS = vtkSmartPointer<vtkMatrix4x4>::New();
  
  // TODO: Can halve the size if there's just a function to handle adding a transform point versus an image point
  // TODO: Vectors can be interpolated instead of being bound to the reference volume resolution
  //Pre-Process Transform Input
  if (inputIsTransform)
  {
    vtkMRMLTransformNode* inputTransformNode = vtkMRMLTransformNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    vtkGeneralTransform* inputTransform = inputTransformNode->GetTransformToParent();
    
    vtkMRMLVolumeNode* referenceVolumeNode = vtkMRMLVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetReferenceVolumeNodeID()));
    referenceVolumeNode->GetIJKToRASMatrix(IJKToRAS);
    
    const int *extent = referenceVolumeNode->GetImageData()->GetExtent();
  
    int iExtent = extent[1]-extent[0];
    int jExtent = extent[3]-extent[2];
    int kExtent = extent[5]-extent[4];
    int numPts = (iExtent+1)*(jExtent+1)*(kExtent+1);

    vtkSmartPointer<vtkMinimalStandardRandomSequence> random = vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();
    random->SetSeed(seed);
    
    int neededValuesAdded = 0;
    int rangeLeft = 0;
    int neededValuesLeft= 0;
    int currentPointID = 0;
    
    double fixedPoint_RAS[4] = {0,0,0,0};
    double movingPoint_RAS[4] = {0,0,0,0};  
    double computedVector_RAS[3] = {0,0,0};
    
    double vMag = 0;
    for(int k = extent[4]; k <= extent[5] && (neededValuesAdded < pointMax); k++)
    {
      for(int j = extent[2]; j <= extent[3] && (neededValuesAdded < pointMax); j++)
      {
        for(int i = extent[0]; i <= extent[1] && (neededValuesAdded < pointMax); i++)
        {
          random->Next();
          neededValuesLeft = pointMax - neededValuesAdded;
          currentPointID = (i+j*(iExtent+1)+k*(iExtent+1)*(jExtent+1));
          rangeLeft = numPts - currentPointID;
          
          if ((random->GetRangeValue(0, rangeLeft)) < neededValuesLeft)
          {
            //Calculate corresponding vector
            fixedPoint_RAS[0]=i; fixedPoint_RAS[1]=j; fixedPoint_RAS[2]=k; fixedPoint_RAS[3]=1; 
            IJKToRAS->MultiplyPoint(fixedPoint_RAS, fixedPoint_RAS);
            inputTransform->TransformPoint(fixedPoint_RAS, movingPoint_RAS);
            computedVector_RAS[0] = movingPoint_RAS[0] - fixedPoint_RAS[0];
            computedVector_RAS[1] = movingPoint_RAS[1] - fixedPoint_RAS[1];
            computedVector_RAS[2] = movingPoint_RAS[2] - fixedPoint_RAS[2];
            vMag = vtkMath::Norm(computedVector_RAS);
            
            //Ignore the point and vector if outside magnitude bounds
            if (vMag >= min && vMag <= max)
            {
              //Add point
              points->InsertNextPoint(fixedPoint_RAS[0], fixedPoint_RAS[1], fixedPoint_RAS[2]);
              
              //Add corresponding vector
              vectors->InsertNextTuple3(movingPoint_RAS[0] - fixedPoint_RAS[0], movingPoint_RAS[1] - fixedPoint_RAS[1], movingPoint_RAS[2] - fixedPoint_RAS[2]);
            }
            neededValuesAdded++;
          }
        }
      }
    }
  }
  //Pre-Process Vector Volume Input
  else
  {
    vtkMRMLVectorVolumeNode* inputVectorVolumeNode = vtkMRMLVectorVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    inputVectorVolumeNode->GetIJKToRASMatrix(IJKToRAS);
    
    vtkImageData* inputDeformationField = inputVectorVolumeNode->GetImageData();
    vtkPointData* deformationFieldPointData = inputDeformationField->GetPointData();
    vtkFloatArray* deformationFieldVectors = vtkFloatArray::SafeDownCast(deformationFieldPointData->GetVectors());
    if (deformationFieldVectors == NULL)
    {
      deformationFieldVectors = vtkFloatArray::SafeDownCast(deformationFieldPointData->GetScalars());
    }
    
    const int *extent = inputDeformationField->GetExtent();
    
    int iExtent = extent[1]-extent[0];
    int jExtent = extent[3]-extent[2];
    int kExtent = extent[5]-extent[4];
    int numPts = (iExtent+1)*(jExtent+1)*(kExtent+1);    
    
    vtkSmartPointer<vtkMinimalStandardRandomSequence> random = vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();
    random->SetSeed(seed);
    
    int neededValuesAdded = 0;
    int rangeLeft = 0;
    int neededValuesLeft= 0;
    int currentPointID = 0;
    
    double fixedPoint_RAS[4] = {0,0,0,0};
    double* chosenVector = NULL;
    
    double vMag = 0;
    for(int k = extent[4]; k <= extent[5] && (neededValuesAdded < pointMax); k++)
    {
      for(int j = extent[2]; j <= extent[3] && (neededValuesAdded < pointMax); j++)
      {
        for(int i = extent[0]; i <= extent[1] && (neededValuesAdded < pointMax); i++)
        {
          random->Next();
          neededValuesLeft = pointMax - neededValuesAdded;
          currentPointID = (i+j*(iExtent+1)+k*(iExtent+1)*(jExtent+1));
          rangeLeft = numPts - currentPointID;
          
          if ((random->GetRangeValue(0, rangeLeft)) < neededValuesLeft)
          {
            //Calculate corresponding vector
            //Just use Applying Matrix thing
            fixedPoint_RAS[0]=i; fixedPoint_RAS[1]=j; fixedPoint_RAS[2]=k; fixedPoint_RAS[3]=1; 
            IJKToRAS->MultiplyPoint(fixedPoint_RAS, fixedPoint_RAS);
            chosenVector = deformationFieldVectors->GetTuple3(currentPointID);
            vMag = vtkMath::Norm(chosenVector, 3);
            
            //Ignore the point and vector if outside magnitude bounds
            if (vMag >= min && vMag <= max)
            {
              //Add point
              points->InsertNextPoint(fixedPoint_RAS[0], fixedPoint_RAS[1], fixedPoint_RAS[2]);
              
              //Add corresponding vector
              vectors->InsertNextTuple3(chosenVector[0], chosenVector[1], chosenVector[2]);
            }
            neededValuesAdded++;
          }
        }
      }
    }
  }

  pointSet->SetPoints(points);
  vtkPointData* pointData = pointSet->GetPointData();
  pointData->SetVectors(vectors);
}

//Grid Visualization
//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::GridVisualization(vtkImageData* field, vtkPolyData* output)
{
  const int subdivision = 1;
  int lineSpacing = this->TransformVisualizerNode->GetGridSpacingMM();
  
  vtkSmartPointer<vtkImageResample> resampled = vtkSmartPointer<vtkImageResample>::New();
  resampled->SetInput(field);
  resampled->SetAxisOutputSpacing(0, lineSpacing / subdivision);
  resampled->SetAxisOutputSpacing(1, lineSpacing / subdivision);
  resampled->SetAxisOutputSpacing(2, lineSpacing / subdivision);
  resampled->Update();

  vtkSmartPointer<vtkImageData> resampledField = vtkSmartPointer<vtkImageData>::New();
  resampledField = resampled->GetOutput();
  
  const double* origin = resampledField->GetOrigin();
  const double* spacing = resampledField->GetSpacing();
  const int* dimensions = resampledField->GetDimensions();  
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(dimensions[0]*dimensions[1]*dimensions[2]);
  vtkSmartPointer<vtkCellArray> grid = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
  
  int i,j,k;
  for (k = 0; k < dimensions[2]; k++)
  {
    for (j = 0; j < dimensions[1]; j++)
    {
      for (i = 0; i < dimensions[0]; i++)
      {
        points->SetPoint(i + j*dimensions[0] + k*dimensions[0]*dimensions[1],
          origin[0] + i + ((spacing[0]-1)*i),origin[1] + j + ((spacing[1]-1)*j),origin[2] + k + ((spacing[2]-1)*k));
      }
    }
  }
  
  for (k = 0; k < dimensions[2]; k+=subdivision)
  {
    for (j = 0; j < dimensions[1]; j+=subdivision)
    {
      for (i = 0; i < ((dimensions[0]-1)-((dimensions[0]-1)%subdivision)); i++)
      {
        line->GetPointIds()->SetId(0, (i) + (j*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
        line->GetPointIds()->SetId(1, (i+1) + (j*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
        grid->InsertNextCell(line);
      }
    }
  }

  for (k = 0; k < dimensions[2]; k+=subdivision)
  {
    for (j = 0; j < ((dimensions[1]-1)-((dimensions[1]-1)%subdivision)); j++)
    {
      for (i = 0; i < dimensions[0]; i+=subdivision)
      {
        line->GetPointIds()->SetId(0, (i) + ((j)*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
        line->GetPointIds()->SetId(1, (i) + ((j+1)*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
        grid->InsertNextCell(line);
      }
    }
  }
  
  for (k = 0; k < ((dimensions[2]-1)-((dimensions[2]-1)%subdivision)); k++)
  {
    for (j = 0; j < dimensions[1]; j+=subdivision)
    {
      for (i = 0; i < dimensions[0]; i+=subdivision)
      {
        line->GetPointIds()->SetId(0, (i) + ((j)*dimensions[0]) + ((k)*dimensions[0]*dimensions[1]));
        line->GetPointIds()->SetId(1, (i) + ((j)*dimensions[0]) + ((k+1)*dimensions[0]*dimensions[1]));
        grid->InsertNextCell(line);
      }
    }
  }  
  
  resampledField->GetPointData()->SetActiveVectors("ImageScalars");
  vtkSmartPointer<vtkFloatArray> vectorArray = vtkFloatArray::SafeDownCast(resampledField->GetPointData()->GetVectors());
  
  vtkSmartPointer<vtkVectorNorm> norm = vtkSmartPointer<vtkVectorNorm>::New();
  norm->SetInputConnection(resampledField->GetProducerPort());
  norm->Update();
  
  vtkSmartPointer<vtkFloatArray> vectorMagnitude = vtkFloatArray::SafeDownCast(norm->GetOutput()->GetPointData()->GetScalars());
  
  vtkSmartPointer<vtkPolyData> polygrid = vtkSmartPointer<vtkPolyData>::New();
  polygrid->SetPoints(points);
  polygrid->SetLines(grid);
  polygrid->GetPointData()->AddArray(vectorArray);
  polygrid->GetPointData()->SetActiveVectors(vectorArray->GetName());
  
  vtkSmartPointer<vtkWarpVector> warp = vtkSmartPointer<vtkWarpVector>::New();
  warp->SetInputConnection(polygrid->GetProducerPort());
  warp->SetScaleFactor(this->TransformVisualizerNode->GetGridScale());
  warp->Update();
  
  output->ShallowCopy(warp->GetPolyDataOutput());
  output->Update();
  output->GetPointData()->AddArray(vectorMagnitude);
  vectorMagnitude->SetName("VectorMagnitude");
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::GridSliceVisualization(vtkImageData* inputField, vtkPolyData* output, vtkSmartPointer<vtkMatrix4x4> rasToIjkDirections)
{
  vtkSmartPointer<vtkImageData> field = vtkSmartPointer<vtkImageData>::New();
  field->DeepCopy(inputField);
  
  vtkSmartPointer<vtkRibbonFilter> ribbon = vtkSmartPointer<vtkRibbonFilter>::New();
  float sliceNormal[3] = {0,0,0};
  double width = 1;

  vtkSmartPointer<vtkMRMLSliceNode> sliceNode = NULL;
  vtkSmartPointer<vtkMRMLNode> node = this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetGridSliceNodeID());
  if (node != NULL)
  {
    sliceNode = vtkMRMLSliceNode::SafeDownCast(node);
  }
  
  vtkSmartPointer<vtkMatrix4x4> sliceToIjk = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMatrix4x4::Multiply4x4(rasToIjkDirections,sliceNode->GetSliceToRAS(),sliceToIjk);
  
  sliceNormal[0] = sliceToIjk->GetElement(0,2);
  sliceNormal[1] = sliceToIjk->GetElement(1,2);
  sliceNormal[2] = sliceToIjk->GetElement(2,2);
  
  const int subdivision = 1;
  int lineSpacing = this->TransformVisualizerNode->GetGridSliceSpacingMM();
  
  vtkSmartPointer<vtkImageResample> resampled = vtkSmartPointer<vtkImageResample>::New();
  resampled->SetInput(field);
  if (sliceNormal[0] < EPSILON)
  {
    resampled->SetAxisOutputSpacing(0, lineSpacing / subdivision);
  }
  if (sliceNormal[1] < EPSILON)
  {
    resampled->SetAxisOutputSpacing(1, lineSpacing / subdivision);
  }
  if (sliceNormal[2] < EPSILON)
  {
    resampled->SetAxisOutputSpacing(2, lineSpacing / subdivision);
  }
  resampled->Update();  
  
  vtkSmartPointer<vtkImageData> resampledField = vtkSmartPointer<vtkImageData>::New();
  resampledField = resampled->GetOutput();  
  
  const double* origin = resampledField->GetOrigin();  
  const double* spacing = resampledField->GetSpacing();
  const int* dimensions = resampledField->GetDimensions();  
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(dimensions[0]*dimensions[1]*dimensions[2]);
  vtkSmartPointer<vtkCellArray> grid = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
  
  //TODO: Need to ensure later that orientation of slice lines up with orientation of volume; no guarantee of that currently
  for (int k = 0; k < dimensions[2]; k++)
  {
    for (int j = 0; j < dimensions[1]; j++)
    {
      for (int i = 0; i < dimensions[0]; i++)
      {
        points->SetPoint(i + j*dimensions[0] + k*dimensions[0]*dimensions[1],
          origin[0] + i + ((spacing[0]-1)*i),origin[1] + j + ((spacing[1]-1)*j),origin[2] + k + ((spacing[2]-1)*k));
      }
    }
  }

  //Reformat not supported
  //TODO: Add support for reformat/oblique slices? 
  if (abs(sliceNormal[0]) < EPSILON)
  {
    for (int k = 0; k < dimensions[2]; k+=subdivision)
    {
      for (int j = 0; j < dimensions[1]; j+=subdivision)
      {
        for (int i = 0; i < ((dimensions[0]-1)-((dimensions[0]-1)%subdivision)); i++)
        {
          line->GetPointIds()->SetId(0, (i) + (j*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
          line->GetPointIds()->SetId(1, (i+1) + (j*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
          grid->InsertNextCell(line);
        }
      }
    }
  }
  else
  {
    width = spacing[0];
  }
  if (abs(sliceNormal[1]) < EPSILON)
  {
    for (int k = 0; k < dimensions[2]; k+=subdivision)
    {
      for (int j = 0; j < ((dimensions[1]-1)-((dimensions[1]-1)%subdivision)); j++)
      {
        for (int i = 0; i < dimensions[0]; i+=subdivision)
        {
          line->GetPointIds()->SetId(0, (i) + ((j)*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
          line->GetPointIds()->SetId(1, (i) + ((j+1)*dimensions[0]) + (k*dimensions[0]*dimensions[1]));
          grid->InsertNextCell(line);
        }
      }
    }
  }
  else
  {
    width = spacing[1];
  }    
  if (abs(sliceNormal[2]) < EPSILON)
  {
    for (int k = 0; k < ((dimensions[2]-1)-((dimensions[2]-1)%subdivision)); k++)
    {
      for (int j = 0; j < dimensions[1]; j+=subdivision)
      {
        for (int i = 0; i < dimensions[0]; i+=subdivision)
        {
          line->GetPointIds()->SetId(0, (i) + ((j)*dimensions[0]) + ((k)*dimensions[0]*dimensions[1]));
          line->GetPointIds()->SetId(1, (i) + ((j)*dimensions[0]) + ((k+1)*dimensions[0]*dimensions[1]));
          grid->InsertNextCell(line);
        }
      }
    }
  }
  else
  {
    width = spacing[2];
  }

  resampledField->GetPointData()->SetActiveVectors("ImageScalars");
  
  vtkSmartPointer<vtkVectorNorm> norm = vtkSmartPointer<vtkVectorNorm>::New();
  norm->SetInputConnection(resampledField->GetProducerPort());
  norm->Update();
  
  vtkSmartPointer<vtkFloatArray> vectorMagnitude = vtkFloatArray::SafeDownCast(norm->GetOutput()->GetPointData()->GetScalars());  
  
  //Projection
  float dot = 0
  ;
  float *ptr = (float *)resampledField->GetScalarPointer();
  for(int i = 0; i < resampledField->GetPointData()->GetScalars()->GetNumberOfTuples()*3; i+=3)
  {
    dot = ptr[i]*sliceNormal[0] + ptr[i+1]*sliceNormal[1] + ptr[i+2]*sliceNormal[2];
    ptr[i] =  ptr[i] - dot*sliceNormal[0];
    ptr[i+1] =  ptr[i+1] - dot*sliceNormal[1];
    ptr[i+2] =  ptr[i+2] - dot*sliceNormal[2];
  }
  ribbon->SetDefaultNormal(sliceNormal[0], sliceNormal[1], sliceNormal[2]);

  vtkSmartPointer<vtkFloatArray> vectorArray = vtkFloatArray::SafeDownCast(resampledField->GetPointData()->GetVectors());

  vtkSmartPointer<vtkPolyData> polygrid = vtkSmartPointer<vtkPolyData>::New();
  polygrid->SetPoints(points);
  polygrid->SetLines(grid);
  polygrid->GetPointData()->AddArray(vectorArray);
  polygrid->GetPointData()->SetActiveVectors(vectorArray->GetName());
  
  vtkSmartPointer<vtkWarpVector> warp = vtkSmartPointer<vtkWarpVector>::New();
  warp->SetInputConnection(polygrid->GetProducerPort());
  warp->SetScaleFactor(this->TransformVisualizerNode->GetGridSliceScale());
  warp->Update();
  
  vtkPolyData* polyoutput = warp->GetPolyDataOutput();
  polyoutput->GetPointData()->AddArray(vectorMagnitude);
  vectorMagnitude->SetName("OriginalVectorMagnitude");
  
  // vtkRibbonFilter
  ribbon->SetInput(polyoutput);
  ribbon->SetWidth(width/2);
  ribbon->SetAngle(90);
  ribbon->UseDefaultNormalOn();
  ribbon->Update();
  
  output->ShallowCopy(ribbon->GetOutput());
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::BlockVisualization(bool inputIsTransform, vtkPolyData* output)
{
  vtkSmartPointer<vtkImageData> field = vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> IJKToRAS = vtkSmartPointer<vtkMatrix4x4>::New();
  if (inputIsTransform)
  {
    vtkMRMLTransformNode* inputTransformNode = vtkMRMLTransformNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    vtkGeneralTransform* inputTransform = inputTransformNode->GetTransformToParent();
    
    vtkMRMLVolumeNode* referenceVolumeNode = vtkMRMLVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetReferenceVolumeNodeID()));
    referenceVolumeNode->GetIJKToRASMatrix(IJKToRAS);
    
    const double *spacing = referenceVolumeNode->GetSpacing();
    
    field->DeepCopy(referenceVolumeNode->GetImageData());
    field->SetScalarTypeToFloat();
    field->SetNumberOfScalarComponents(3);
    field->AllocateScalars();
    field->SetSpacing(spacing[0], spacing[1], spacing[2]);
    field->GetPointData()->SetActiveVectors("ImageScalars");
    field->Update();    
    
    const int *extent = referenceVolumeNode->GetImageData()->GetExtent();
    
    double fixedPoint_RAS[4] = {0,0,0,0};
    double movingPoint_RAS[4] = {0,0,0,0};  

    float *ptr = (float *)field->GetScalarPointer();
    for(int k = extent[4]; k < extent[5]+1; k++)
    {
      for(int j = extent[2]; j < extent[3]+1; j++)
      {
        for(int i = extent[0]; i < extent[1]+1; i++)
        {
          fixedPoint_RAS[0]=i; fixedPoint_RAS[1]=j; fixedPoint_RAS[2]=k; fixedPoint_RAS[3]=1; 
          IJKToRAS->MultiplyPoint(fixedPoint_RAS, fixedPoint_RAS);
          inputTransform->TransformPoint(fixedPoint_RAS, movingPoint_RAS);
          ptr[(i + j*(extent[1]+1) + k*(extent[1]+1)*(extent[3]+1))*3] = movingPoint_RAS[0] - fixedPoint_RAS[0];
          ptr[(i + j*(extent[1]+1) + k*(extent[1]+1)*(extent[3]+1))*3 + 1] = movingPoint_RAS[1] - fixedPoint_RAS[1];
          ptr[(i + j*(extent[1]+1) + k*(extent[1]+1)*(extent[3]+1))*3 + 2] = movingPoint_RAS[2] - fixedPoint_RAS[2];
        }
      }
    }    
  }
  else
  {
    vtkMRMLVectorVolumeNode* inputVectorVolumeNode = vtkMRMLVectorVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image->DeepCopy(inputVectorVolumeNode->GetImageData());
    
    inputVectorVolumeNode->GetIJKToRASMatrix(IJKToRAS);
    
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();    
    transform->SetMatrix(IJKToRAS);
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetTransform(transform);
    transformFilter->SetInputConnection(image->GetProducerPort());
    
    field->ShallowCopy(transformFilter->GetOutput());
  }

  vtkSmartPointer<vtkVectorNorm> norm = vtkSmartPointer<vtkVectorNorm>::New();
  norm->SetInputConnection(field->GetProducerPort());
  norm->Update();
  vtkSmartPointer<vtkFloatArray> vectorMagnitude = vtkFloatArray::SafeDownCast(norm->GetOutput()->GetPointData()->GetScalars());

  vtkSmartPointer<vtkWarpVector> warp = vtkSmartPointer<vtkWarpVector>::New();
  warp->SetInputConnection(field->GetProducerPort());
  warp->SetScaleFactor(this->TransformVisualizerNode->GetBlockScale());

  //TODO: Current method of generating polydata is very inefficient but avoids bugs with possible extreme cases. Better method to be implemented.
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputConnection(warp->GetOutputPort());
  
  vtkPolyData* polyoutput = geometryFilter->GetOutput();
  polyoutput->Update();
  polyoutput->GetPointData()->AddArray(vectorMagnitude);
  vectorMagnitude->SetName("VectorMagnitude");
  
  if (this->TransformVisualizerNode->GetBlockDisplacementCheck())
  {
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInput(polyoutput);  
    normals->Update();
    normals->GetOutput()->GetPointData()->SetVectors(field->GetPointData()->GetScalars());
    
    vtkSmartPointer<vtkVectorDot> dot = vtkSmartPointer<vtkVectorDot>::New();
    dot->SetInput(normals->GetOutput());
    dot->Update();
    vtkSmartPointer<vtkFloatArray> vectorDot = vtkFloatArray::SafeDownCast(dot->GetOutput()->GetPointData()->GetScalars());
    
    normals->GetOutput()->GetPointData()->AddArray(vectorDot);
    vectorDot->SetName("VectorDot");
    
    output->ShallowCopy(normals->GetOutput());
  }

  output->ShallowCopy(polyoutput);
}

//----------------------------------------------------------------------------
void vtkSlicerTransformVisualizerLogic::ContourVisualization(bool inputIsTransform, vtkPolyData* output)
{
  vtkSmartPointer<vtkImageData> scalarMagnitudes = vtkSmartPointer<vtkImageData>::New();
  if (inputIsTransform)
  {
  }
  else
  {
    vtkMRMLVectorVolumeNode* inputVectorVolumeNode = vtkMRMLVectorVolumeNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(this->TransformVisualizerNode->GetInputNodeID()));
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image->DeepCopy(inputVectorVolumeNode->GetImageData());
    
    vtkSmartPointer<vtkVectorNorm> norm = vtkSmartPointer<vtkVectorNorm>::New();
    norm->SetInputConnection(image->GetProducerPort());
    norm->Update();
    scalarMagnitudes->DeepCopy(norm->GetOutput());
  }
   
  vtkSmartPointer<vtkMarchingCubes> iso = vtkSmartPointer<vtkMarchingCubes>::New();
	iso->SetInputConnection(scalarMagnitudes->GetProducerPort());
  iso->ComputeScalarsOn();
  iso->ComputeNormalsOff();
  iso->ComputeGradientsOff();
  double* values = this->TransformVisualizerNode->GetContourValues();
  iso->SetNumberOfContours(this->TransformVisualizerNode->GetContourNumber());
  for(int i = 0; i < this->TransformVisualizerNode->GetContourNumber(); i++)
  {
    iso->SetValue(i, values[i]);
  }
  iso->Update();
  
  vtkSmartPointer<vtkDecimatePro> decimator = vtkSmartPointer<vtkDecimatePro>::New();
  decimator->SetInputConnection(iso->GetOutputPort());
  decimator->SetFeatureAngle(60);
  decimator->SplittingOff();
  decimator->PreserveTopologyOn();
  decimator->SetMaximumError(1);
  decimator->SetTargetReduction(this->TransformVisualizerNode->GetContourDecimation());
  decimator->Update();

  decimator->GetOutput()->GetPointData()->GetScalars()->SetName("VectorMagnitude");
  
  output->ShallowCopy(decimator->GetOutput());
}
