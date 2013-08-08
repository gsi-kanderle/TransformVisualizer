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

// Qt includes
#include <QDebug>
#include <QtCore>
#include <QtGui>

// SlicerQt includes
#include "qSlicerTransformVisualizerModuleWidget.h"
#include "ui_qSlicerTransformVisualizerModule.h"
#include <qSlicerApplication.h>

// DeformationVisualizer includes
#include "vtkSlicerTransformVisualizerLogic.h"
#include "vtkMRMLTransformVisualizerNode.h"

// MMRL includes
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLLinearTransformNode.h>
#include <vtkMRMLBSplineTransformNode.h>
#include <vtkMRMLGridTransformNode.h>
#include <vtkMRMLSliceNode.h>

// VTK includes
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkGeneralTransform.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_TransformVisualizer
// TODO: Keeping private for now until after fixes and enhancements
class qSlicerTransformVisualizerModuleWidgetPrivate: public Ui_qSlicerTransformVisualizerModule{
  Q_DECLARE_PUBLIC(qSlicerTransformVisualizerModuleWidget);
protected:
  qSlicerTransformVisualizerModuleWidget* const q_ptr;
public:
  qSlicerTransformVisualizerModuleWidgetPrivate(qSlicerTransformVisualizerModuleWidget& object);
  vtkSlicerTransformVisualizerLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerTransformVisualizerModuleWidgetPrivate methods
//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModuleWidgetPrivate::qSlicerTransformVisualizerModuleWidgetPrivate(qSlicerTransformVisualizerModuleWidget& object) : q_ptr(&object)
{
}

vtkSlicerTransformVisualizerLogic* qSlicerTransformVisualizerModuleWidgetPrivate::logic() const {
  Q_Q( const qSlicerTransformVisualizerModuleWidget );
  return vtkSlicerTransformVisualizerLogic::SafeDownCast( q->logic() );
}


//-----------------------------------------------------------------------------
// qSlicerTransformVisualizerModuleWidget methods
//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModuleWidget::qSlicerTransformVisualizerModuleWidget(QWidget* _parent) : Superclass( _parent ) , d_ptr(new qSlicerTransformVisualizerModuleWidgetPrivate(*this))
{
}

//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModuleWidget::~qSlicerTransformVisualizerModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setMRMLScene(vtkMRMLScene* scene)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  this->Superclass::setMRMLScene(scene);

  // Find parameters node or create it if there is no one in the scene
  if (scene &&  d->logic()->GetTransformVisualizerNode() == 0)
  {
    vtkMRMLNode* node = scene->GetNthNodeByClass(0, "vtkMRMLTransformVisualizerNode");
    if (node){
      this->setTransformVisualizerParametersNode(vtkMRMLTransformVisualizerNode::SafeDownCast(node));
    }
  }
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::onSceneImportedEvent()
{
  this->onEnter();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::enter()
{
  this->onEnter();
  this->Superclass::enter();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::onEnter()
{
  if (!this->mrmlScene())
  {
    return;
  }

  Q_D(qSlicerTransformVisualizerModuleWidget);

  if (d->logic() == NULL)
  {
    return;
  }

  //Check for existing parameter node
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();

  if (pNode == NULL)
  {
    vtkMRMLNode* node = this->mrmlScene()->GetNthNodeByClass(0, "vtkMRMLTransformVisualizerNode");
    if (node)
    {
      pNode = vtkMRMLTransformVisualizerNode::SafeDownCast(node);
      d->logic()->SetAndObserveTransformVisualizerNode(pNode);
      return;
    }
    else
    {
      vtkSmartPointer<vtkMRMLTransformVisualizerNode> newNode = vtkSmartPointer<vtkMRMLTransformVisualizerNode>::New();
      this->mrmlScene()->AddNode(newNode);
      d->logic()->SetAndObserveTransformVisualizerNode(newNode);
    }
  }
  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setTransformVisualizerParametersNode(vtkMRMLNode *node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = vtkMRMLTransformVisualizerNode::SafeDownCast(node);

  qvtkReconnect( d->logic()->GetTransformVisualizerNode(), pNode, vtkCommand::ModifiedEvent, this, SLOT(updateWidgetFromMRML()));

  d->logic()->SetAndObserveTransformVisualizerNode(pNode);

  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::updateWidgetFromMRML()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();

  if (pNode && this->mrmlScene())
  {
    d->ParameterComboBox->setCurrentNode(pNode);

    if (pNode->GetInputVolumeNodeID())
    {
      d->InputFieldComboBox->setCurrentNodeID(pNode->GetInputVolumeNodeID());
    }
    else
    {
      this->inputVolumeChanged(d->InputFieldComboBox->currentNode());
    }

    if (pNode->GetReferenceVolumeNodeID())
    {
      d->InputReferenceComboBox->setCurrentNodeID(pNode->GetReferenceVolumeNodeID());
    }
    else
    {
      this->referenceVolumeChanged(d->InputReferenceComboBox->currentNode());
    }    

    if (pNode->GetOutputModelNodeID())
    {
      d->OutputModelComboBox->setCurrentNodeID(pNode->GetOutputModelNodeID());
    }
    else
    {
      this->outputModelChanged(d->OutputModelComboBox->currentNode());
    }  

    //Update Visualization Parameters

    // Glyph Parameters
    d->InputGlyphPointMax->setValue(pNode->GetGlyphPointMax());
    d->InputGlyphSeed->setValue(pNode->GetGlyphSeed());
    d->InputGlyphScale->setValue(pNode->GetGlyphScale());
    d->InputGlyphScaleDirectional->setChecked(pNode->GetGlyphScaleDirectional());
    d->InputGlyphScaleIsotropic->setChecked(pNode->GetGlyphScaleIsotropic());
    d->InputGlyphThreshold->setMaximumValue(pNode->GetGlyphThresholdMax());
    d->InputGlyphThreshold->setMinimumValue(pNode->GetGlyphThresholdMin());
    d->GlyphSourceComboBox->setCurrentIndex(pNode->GetGlyphSourceOption());
    // Arrow Parameters
    d->InputGlyphArrowTipLength->setValue(pNode->GetGlyphArrowTipLength());
    d->InputGlyphArrowTipRadius->setValue(pNode->GetGlyphArrowTipRadius());
    d->InputGlyphArrowShaftRadius->setValue(pNode->GetGlyphArrowShaftRadius());
    d->InputGlyphArrowResolution->setValue(pNode->GetGlyphArrowResolution());
    // Cone Parameters
    d->InputGlyphConeHeight->setValue(pNode->GetGlyphConeHeight());
    d->InputGlyphConeRadius->setValue(pNode->GetGlyphConeRadius());
    d->InputGlyphConeResolution->setValue(pNode->GetGlyphConeResolution());
    // Sphere Parameters
    d->InputGlyphSphereResolution->setValue(pNode->GetGlyphSphereResolution());

    // Grid Parameters
    d->InputGridScale->setValue(pNode->GetGridScale());
    d->InputGridSpacing->setValue(pNode->GetGridSpacingMM());

    // Block Parameters
    d->InputBlockScale->setValue(pNode->GetBlockScale());
    d->InputBlockDisplacementCheck->setChecked(pNode->GetBlockDisplacementCheck());

    // Contour Parameters
    d->InputContourNumber->setValue(pNode->GetContourNumber());
    d->InputContourRange->setMaximumValue(pNode->GetContourMax());
    d->InputContourRange->setMinimumValue(pNode->GetContourMin());
    d->InputContourDecimation->setValue(pNode->GetContourDecimation());

    // Glyph Slice Parameters
    if (pNode->GetGlyphSliceNodeID())
    {
      d->GlyphSliceComboBox->setCurrentNodeID(pNode->GetGlyphSliceNodeID());
    }
    else
    {
      this->setGlyphSliceNode(d->GlyphSliceComboBox->currentNode());
    }  
    d->InputGlyphSlicePointMax->setValue(pNode->GetGlyphSlicePointMax());
    d->InputGlyphSliceThreshold->setMaximumValue(pNode->GetGlyphSliceThresholdMax());
    d->InputGlyphSliceThreshold->setMinimumValue(pNode->GetGlyphSliceThresholdMin());
    d->InputGlyphSliceScale->setValue(pNode->GetGlyphSliceScale());
    d->InputGlyphSliceSeed->setValue(pNode->GetGlyphSliceSeed());

    // Grid Slice Parameters
    if (pNode->GetGridSliceNodeID())
    {
      d->GridSliceComboBox->setCurrentNodeID(pNode->GetGridSliceNodeID());
    }
    else
    {
      this->setGridSliceNode(d->GridSliceComboBox->currentNode());
    }  
    d->InputGridSliceScale->setValue(pNode->GetGridSliceScale());
    d->InputGridSliceSpacing->setValue(pNode->GetGridSliceSpacingMM());
  }
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::onLogicModified()
{
  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::inputVolumeChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  // TODO: Move into updatefrommrml?
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene() || !node)
  {
    d->ApplyButton->setEnabled(false);
    d->VolumeDisabledLabel->show();
    return;
  }

  d->ApplyButton->setEnabled(true);
  d->VolumeDisabledLabel->hide();  

  pNode->DisableModifiedEventOn();
  pNode->SetAndObserveInputVolumeNodeID(node->GetID());
  pNode->DisableModifiedEventOff();
  
  double* range;
  
  // What to do if there is more than one array? Would there be more than one array?
  if (strcmp(node->GetClassName(), "vtkMRMLVectorVolumeNode") == 0)
  {
    d->InputReferenceComboBox->setEnabled(false);
    range = vtkMRMLVectorVolumeNode::SafeDownCast(node)->GetImageData()->GetPointData()->GetScalars()->GetRange(-1);
  }
  else if (strcmp(node->GetClassName(), "vtkMRMLLinearTransformNode") == 0 || 
    strcmp(node->GetClassName(), "vtkMRMLBSplineTransformNode") == 0 ||
    strcmp(node->GetClassName(), "vtkMRMLGridTransformNode") == 0)
  {
    d->InputReferenceComboBox->setEnabled(true);
    vtkSmartPointer<vtkMRMLVolumeNode> referenceVolumeNode = vtkMRMLVolumeNode::SafeDownCast(this->mrmlScene()->GetNodeByID(pNode->GetReferenceVolumeNodeID()));

    if (referenceVolumeNode == NULL)
    {
      return;
    }

    //TODO: Remake progress dialog and add detail (update progress from actual steps occurring in logic)
    QProgressDialog *convertProgress =  new QProgressDialog(qSlicerApplication::application()->mainWindow());
    convertProgress->setCancelButton(0);
    convertProgress->setModal(true);
    convertProgress->setMinimumDuration(100);
    convertProgress->show();
    convertProgress->setLabelText("Converting transform to vector volume...");
    
    convertProgress->setValue(20);
    d->logic()->GenerateDeformationField();
    
    convertProgress->setValue(80);
    range = d->logic()->GetFieldRange();
    
    convertProgress->setValue(100);
    delete convertProgress;
  }
  else
  {
    return;
  }

  pNode->SetGlyphThresholdMin(range[0]);
  d->InputGlyphThreshold->setMinimum(range[0]);
  d->InputGlyphThreshold->setMinimumValue(range[0]);
  pNode->SetGlyphThresholdMax(range[1]);
  d->InputGlyphThreshold->setMaximum(range[1]);
  d->InputGlyphThreshold->setMaximumValue(range[1]);

  pNode->SetContourMin(range[0]);
  d->InputContourRange->setMinimum(range[0]);
  d->InputContourRange->setMinimumValue(range[0]);
  pNode->SetContourMax(range[1]);
  d->InputContourRange->setMaximum(range[1]);
  d->InputContourRange->setMaximumValue(range[1]);

  pNode->SetGlyphSliceThresholdMin(range[0]);
  d->InputGlyphSliceThreshold->setMinimum(range[0]);
  d->InputGlyphSliceThreshold->setMinimumValue(range[0]);
  pNode->SetGlyphSliceThresholdMax(range[1]);
  d->InputGlyphSliceThreshold->setMaximum(range[1]);
  d->InputGlyphSliceThreshold->setMaximumValue(range[1]);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::referenceVolumeChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene() || !node)
  {
    return;
  }

  pNode->DisableModifiedEventOn();
  pNode->SetAndObserveReferenceVolumeNodeID(node->GetID());
  pNode->DisableModifiedEventOff();
  
  vtkSmartPointer<vtkMRMLTransformNode> inputVolumeNode = vtkMRMLTransformNode::SafeDownCast(this->mrmlScene()->GetNodeByID(pNode->GetInputVolumeNodeID()));
  if (inputVolumeNode == NULL)
  {
    return;
  }
  
  //TODO: Remake progress dialog and add detail (update progress from actual steps occurring in logic)
    QProgressDialog *convertProgress =  new QProgressDialog(qSlicerApplication::application()->mainWindow());
    convertProgress->setCancelButton(0);
    convertProgress->setModal(true);
    convertProgress->setMinimumDuration(100);
    convertProgress->show();
  convertProgress->setLabelText("Converting transform to vector volume...");
  
  convertProgress->setValue(20);
  d->logic()->GenerateDeformationField();
  
  convertProgress->setValue(80);
  double* range = d->logic()->GetFieldRange();
  
  convertProgress->setValue(100);
  delete convertProgress;

  pNode->SetGlyphThresholdMin(range[0]);
  d->InputGlyphThreshold->setMinimum(range[0]);
  d->InputGlyphThreshold->setMinimumValue(range[0]);
  pNode->SetGlyphThresholdMax(range[1]);
  d->InputGlyphThreshold->setMaximum(range[1]);
  d->InputGlyphThreshold->setMaximumValue(range[1]);

  pNode->SetContourMin(range[0]);
  d->InputContourRange->setMinimum(range[0]);
  d->InputContourRange->setMinimumValue(range[0]);
  pNode->SetContourMax(range[1]);
  d->InputContourRange->setMaximum(range[1]);
  d->InputContourRange->setMaximumValue(range[1]);

  pNode->SetGlyphSliceThresholdMin(range[0]);
  d->InputGlyphSliceThreshold->setMinimum(range[0]);
  d->InputGlyphSliceThreshold->setMinimumValue(range[0]);
  pNode->SetGlyphSliceThresholdMax(range[1]);
  d->InputGlyphSliceThreshold->setMaximum(range[1]);
  d->InputGlyphSliceThreshold->setMaximumValue(range[1]);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::outputModelChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene() || !node)
  {
    d->ApplyButton->setEnabled(false);
    d->ModelDisabledLabel->show();
    return;
  }

  d->ApplyButton->setEnabled(true);
  d->ModelDisabledLabel->hide();

  pNode->DisableModifiedEventOn();
  pNode->SetAndObserveOutputModelNodeID(node->GetID());
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::updateSourceOptions(int option)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();

  if (option == d->logic()->ARROW_3D)
  {
    d->ArrowSourceOptions->setEnabled(true);
    d->ArrowSourceOptions->setVisible(true);
    d->ConeSourceOptions->setEnabled(false);
    d->ConeSourceOptions->setVisible(false);
    d->SphereSourceOptions->setEnabled(false);
    d->SphereSourceOptions->setVisible(false);

    if (!pNode || !this->mrmlScene())
  {
    return;
  }
    pNode->DisableModifiedEventOn();
    pNode->SetGlyphScaleDirectional(true);
    pNode->SetGlyphScaleIsotropic(false);
    pNode->DisableModifiedEventOff();
  }
  else if (option == d->logic()->CONE_3D)
  {
    d->ArrowSourceOptions->setEnabled(false);
    d->ArrowSourceOptions->setVisible(false);
    d->ConeSourceOptions->setEnabled(true);
    d->ConeSourceOptions->setVisible(true);
    d->SphereSourceOptions->setEnabled(false);
    d->SphereSourceOptions->setVisible(false);  

    if (!pNode || !this->mrmlScene())
  {
    return;
  }
    pNode->DisableModifiedEventOn();
    pNode->SetGlyphScaleDirectional(true);
    pNode->SetGlyphScaleIsotropic(false);
    pNode->DisableModifiedEventOff();
  }
  else if (option == d->logic()->SPHERE_3D)
  {
    d->ArrowSourceOptions->setEnabled(false);
    d->ArrowSourceOptions->setVisible(false);
    d->ConeSourceOptions->setEnabled(false);
    d->ConeSourceOptions->setVisible(false);
    d->SphereSourceOptions->setEnabled(true);
    d->SphereSourceOptions->setVisible(true);

    if (!pNode || !this->mrmlScene())
  {
    return;
  }
    pNode->DisableModifiedEventOn();
    pNode->SetGlyphScaleDirectional(false);
    pNode->SetGlyphScaleIsotropic(true);
    pNode->DisableModifiedEventOff();
  }

  d->InputGlyphScaleDirectional->setChecked(pNode->GetGlyphScaleDirectional());
  d->InputGlyphScaleIsotropic->setChecked(pNode->GetGlyphScaleIsotropic());
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::visualize()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  if (d->InputFieldComboBox->currentNodeID() != NULL && d->OutputModelComboBox->currentNodeID() != NULL)
  {
  //TODO: Remake progress dialog and add detail (update progress from actual steps occurring in logic)
    QProgressDialog *visualizeProgress =  new QProgressDialog(qSlicerApplication::application()->mainWindow());
    visualizeProgress->setCancelButton(0);
    visualizeProgress->setModal(true);
    visualizeProgress->setMinimumDuration(100); //will matter a bit more after progress dialog is remade
    visualizeProgress->show();
  visualizeProgress->setLabelText("Processing...");
    visualizeProgress->setValue(0);
  
  if (d->GlyphToggle->isChecked())
  {
    visualizeProgress->setLabelText("Creating glyphs...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_GLYPH_3D);
    }
    else if (d->GridToggle->isChecked()){
    visualizeProgress->setLabelText("Creating grid...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_GRID_3D);
    }
    else if (d->BlockToggle->isChecked()){
    visualizeProgress->setLabelText("Creating block...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_BLOCK_3D);
    }
    else if (d->ContourToggle->isChecked()){
    visualizeProgress->setLabelText("Creating contours...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_CONTOUR_3D);
    }    
    else if (d->GlyphSliceToggle->isChecked()){
    visualizeProgress->setLabelText("Creating glyphs for slice view...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_GLYPH_2D);
    }
    else if (d->GridSliceToggle->isChecked()){
    visualizeProgress->setLabelText("Creating grid for slice view...");
    visualizeProgress->setValue(20);
    d->logic()->CreateVisualization(d->logic()->VIS_MODE_GRID_2D);
    }
  visualizeProgress->setValue(100);
  delete visualizeProgress;
  }
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setup()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  connect(d->ParameterComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(setTransformVisualizerParametersNode(vtkMRMLNode*)));

  connect(d->InputFieldComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(inputVolumeChanged(vtkMRMLNode*)));
  connect(d->InputReferenceComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(referenceVolumeChanged(vtkMRMLNode*)));
  connect(d->OutputModelComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(outputModelChanged(vtkMRMLNode*)));

  // Glyph Parameters
  connect(d->InputGlyphPointMax, SIGNAL(valueChanged(double)), this, SLOT(setGlyphPointMax(double)));
  connect(d->InputGlyphScale, SIGNAL(valueChanged(double)), this, SLOT(setGlyphScale(double)));
  connect(d->InputGlyphScaleDirectional, SIGNAL(toggled(bool)), this, SLOT(setGlyphScaleDirectional(bool)));
  connect(d->InputGlyphScaleIsotropic, SIGNAL(toggled(bool)), this, SLOT(setGlyphScaleIsotropic(bool)));
  connect(d->InputGlyphThreshold, SIGNAL(valuesChanged(double, double)), this, SLOT(setGlyphThreshold(double, double)));
  connect(d->GenerateSeedButton, SIGNAL(clicked()), this, SLOT(setSeed()));
  connect(d->InputGlyphSeed, SIGNAL(valueChanged(int)), this, SLOT(setGlyphSeed(int)));  
  connect(d->GlyphSourceComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setGlyphSourceOption(int)));
  // Arrow Parameters
  connect(d->InputGlyphArrowTipLength, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowTipLength(double)));
  connect(d->InputGlyphArrowTipRadius, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowTipRadius(double)));
  connect(d->InputGlyphArrowShaftRadius, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowShaftRadius(double)));  
  connect(d->InputGlyphArrowResolution, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowResolution(double)));
  // Cone Parameters
  connect(d->InputGlyphConeHeight, SIGNAL(valueChanged(double)), this, SLOT(setGlyphConeHeight(double)));
  connect(d->InputGlyphConeRadius, SIGNAL(valueChanged(double)), this, SLOT(setGlyphConeRadius(double)));
  connect(d->InputGlyphConeResolution, SIGNAL(valueChanged(double)), this, SLOT(setGlyphConeResolution(double)));
  // Sphere Parameters
  connect(d->InputGlyphSphereResolution, SIGNAL(valueChanged(double)), this, SLOT(setGlyphSphereResolution(double)));

  // Grid Parameters
  connect(d->InputGridScale, SIGNAL(valueChanged(double)), this, SLOT(setGridScale(double)));
  connect(d->InputGridSpacing, SIGNAL(valueChanged(double)), this, SLOT(setGridSpacingMM(double)));

  // Block Parameters  
  connect(d->InputBlockScale, SIGNAL(valueChanged(double)), this, SLOT(setBlockScale(double)));
  connect(d->InputBlockDisplacementCheck, SIGNAL(stateChanged(int)), this, SLOT(setBlockDisplacementCheck(int)));

  // Contour Parameters
  connect(d->InputContourNumber, SIGNAL(valueChanged(double)), this, SLOT(setContourNumber(double)));
  connect(d->InputContourRange, SIGNAL(valuesChanged(double, double)), this, SLOT(setContourRange(double, double)));
  connect(d->InputContourDecimation, SIGNAL(valueChanged(double)), this, SLOT(setContourDecimation(double)));

  // Glyph Slice Parameters
  connect(d->GlyphSliceComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(setGlyphSliceNode(vtkMRMLNode*)));
  connect(d->InputGlyphSlicePointMax, SIGNAL(valueChanged(double)), this, SLOT(setGlyphSlicePointMax(double)));  
  connect(d->InputGlyphSliceThreshold, SIGNAL(valuesChanged(double, double)), this, SLOT(setGlyphSliceThreshold(double, double)));
  connect(d->InputGlyphSliceScale, SIGNAL(valueChanged(double)), this, SLOT(setGlyphSliceScale(double)));
  connect(d->InputGlyphSliceSeed, SIGNAL(valueChanged(int)), this, SLOT(setGlyphSliceSeed(int)));  
  connect(d->GenerateSeedButton2, SIGNAL(clicked()), this, SLOT(setSeed2()));

  // Grid Slice Parameters
  connect(d->GridSliceComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(setGridSliceNode(vtkMRMLNode*)));
  connect(d->InputGridSliceScale, SIGNAL(valueChanged(double)), this, SLOT(setGridSliceScale(double)));
  connect(d->InputGridSliceSpacing, SIGNAL(valueChanged(double)), this, SLOT(setGridSliceSpacingMM(double)));

  connect(d->ApplyButton, SIGNAL(clicked()), this, SLOT(visualize()));
}

//-----------------------------------------------------------------------------
// Glyph parameters
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphPointMax(double pointMax)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphPointMax(pointMax);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setSeed()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSeed(rand());
  d->InputGlyphSeed->setValue(pNode->GetGlyphSeed());
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSeed(int seed)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSeed(seed);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphScale(double scale)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphScale(scale);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphScaleDirectional(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphScaleDirectional(state);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphScaleIsotropic(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphScaleIsotropic(state);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphThreshold(double min, double max)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphThresholdMin(min);
  pNode->SetGlyphThresholdMax(max);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSourceOption(int option)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSourceOption(option);
  pNode->DisableModifiedEventOff();
  this->updateSourceOptions(option);
}

//-----------------------------------------------------------------------------
// Arrow parameters
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowTipLength(double length)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphArrowTipLength(length);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowTipRadius(double radius)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphArrowTipRadius(radius);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowShaftRadius(double radius)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphArrowShaftRadius(radius);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowResolution(double resolution)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphArrowResolution(resolution);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
// Cone Parameters
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeHeight(double height)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphConeHeight(height);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeRadius(double radius)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphConeRadius(radius);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeResolution(double resolution)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphConeResolution(resolution);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
// Sphere Parameters
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSphereResolution(double resolution)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSphereResolution(resolution);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGridScale(double scale)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGridScale(scale);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGridSpacingMM(double spacing)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGridSpacingMM(spacing);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setBlockScale(double scale)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetBlockScale(scale);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setBlockDisplacementCheck(int state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetBlockDisplacementCheck(state);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setContourNumber(double number)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetContourNumber(number);
  pNode->DisableModifiedEventOff();
}
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setContourRange(double min, double max)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetContourMin(min);
  pNode->SetContourMax(max);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setContourDecimation(double reduction)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetContourDecimation(reduction);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSliceNode(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!node || !pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetAndObserveGlyphSliceNodeID(node->GetID());
  pNode->DisableModifiedEventOff();  
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSlicePointMax(double pointMax)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSlicePointMax(pointMax);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSliceThreshold(double min, double max)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSliceThresholdMin(min);
  pNode->SetGlyphSliceThresholdMax(max);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSliceScale(double scale)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSliceScale(scale);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphSliceSeed(int seed)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSliceSeed(seed);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setSeed2()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGlyphSliceSeed(rand());
  d->InputGlyphSliceSeed->setValue(pNode->GetGlyphSliceSeed());
  pNode->DisableModifiedEventOff();
}


//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGridSliceNode(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!node || !pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetAndObserveGridSliceNodeID(node->GetID());
  pNode->DisableModifiedEventOff();  
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGridSliceScale(double scale)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGridSliceScale(scale);
  pNode->DisableModifiedEventOff();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGridSliceSpacingMM(double spacing)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->DisableModifiedEventOn();
  pNode->SetGridSliceSpacingMM(spacing);
  pNode->DisableModifiedEventOff();
}
