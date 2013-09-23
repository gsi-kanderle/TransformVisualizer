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
  Q_D(qSlicerTransformVisualizerModuleWidget);
  
  if (!this->mrmlScene() || d->logic() == NULL)
  {
    std::cerr << "Error: Unable to initialize module" << std::endl;
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
    }
    else
    {
      vtkSmartPointer<vtkMRMLTransformVisualizerNode> newNode = vtkSmartPointer<vtkMRMLTransformVisualizerNode>::New();
      this->mrmlScene()->AddNode(newNode);
      d->logic()->SetAndObserveTransformVisualizerNode(newNode);
    }
  }
  this->update();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setTransformVisualizerParametersNode(vtkMRMLNode *node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  
  vtkMRMLTransformVisualizerNode* pNode = vtkMRMLTransformVisualizerNode::SafeDownCast(node);
  qvtkReconnect( d->logic()->GetTransformVisualizerNode(), pNode, vtkCommand::ModifiedEvent, this, SLOT(update()));
  d->logic()->SetAndObserveTransformVisualizerNode(pNode);
  this->update();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::update()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (pNode == NULL || this->mrmlScene() == NULL)
  {
    std::cerr << "Error: Unable to update widget" << std::endl;
    return;
  }
  
  // Update widget from MRML
  d->ParameterComboBox->setCurrentNode(pNode);
  d->InputFieldComboBox->setCurrentNodeID(pNode->GetInputNodeID());
  d->InputReferenceComboBox->setCurrentNodeID(pNode->GetReferenceVolumeNodeID());
  d->OutputModelComboBox->setCurrentNodeID(pNode->GetOutputModelNodeID());

  // Update Visualization Parameters
  // Glyph Parameters
  d->InputGlyphPointMax->setValue(pNode->GetGlyphPointMax());
  d->InputGlyphSeed->setValue(pNode->GetGlyphSeed());
  d->InputGlyphScale->setValue(pNode->GetGlyphScale());
  d->InputGlyphThreshold->setMaximumValue(pNode->GetGlyphThresholdMax());
  d->InputGlyphThreshold->setMinimumValue(pNode->GetGlyphThresholdMin());
  d->GlyphSourceComboBox->setCurrentIndex(pNode->GetGlyphSourceOption());
  // Arrow Parameters
  d->InputGlyphArrowScaleDirectional->setChecked(pNode->GetGlyphArrowScaleDirectional());
  d->InputGlyphArrowScaleIsotropic->setChecked(pNode->GetGlyphArrowScaleIsotropic());
  d->InputGlyphArrowTipLength->setValue(pNode->GetGlyphArrowTipLength());
  d->InputGlyphArrowTipRadius->setValue(pNode->GetGlyphArrowTipRadius());
  d->InputGlyphArrowShaftRadius->setValue(pNode->GetGlyphArrowShaftRadius());
  d->InputGlyphArrowResolution->setValue(pNode->GetGlyphArrowResolution());
  // Cone Parameters
  d->InputGlyphConeScaleDirectional->setChecked(pNode->GetGlyphConeScaleDirectional());
  d->InputGlyphConeScaleIsotropic->setChecked(pNode->GetGlyphConeScaleIsotropic());
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
  //d->InputContourNumber->setValue(pNode->GetContourNumber());
  //d->InputContourRange->setMaximumValue(pNode->GetContourMax());
  //d->InputContourRange->setMinimumValue(pNode->GetContourMin());
  d->InputContourDecimation->setValue(pNode->GetContourDecimation());

  // Glyph Slice Parameters
  // Set default to a slice node that exists in the scene
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
 
  this->updateLabels();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::updateLabels()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);

  if (d->InputFieldComboBox->currentNode() == NULL)
  {
    d->ApplyButton->setEnabled(false);
    d->VolumeDisabledLabel->show();
  }
  else
  {
    d->ApplyButton->setEnabled(true);
    d->VolumeDisabledLabel->hide();
    
    if (!d->InputFieldComboBox->currentNode()->IsA("vtkMRMLVectorVolumeNode") && d->InputReferenceComboBox->currentNode() == NULL)
    {
      d->ApplyButton->setEnabled(false);
      d->ReferenceDisabledLabel->show();
    }
    else
    {
      d->ApplyButton->setEnabled(true);
      d->ReferenceDisabledLabel->hide();
    }
  }
  
  if (d->OutputModelComboBox->currentNode() == NULL)
  {
    d->ApplyButton->setEnabled(false);
    d->ModelDisabledLabel->show();
  }
  else
  {
    d->ApplyButton->setEnabled(true);
    d->ModelDisabledLabel->hide();
  }
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::updateGlyphSourceOptions(int sourceOption)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  
  if (sourceOption == d->logic()->ARROW_3D)
  {
    d->ArrowSourceOptions->setEnabled(true);
    d->ArrowSourceOptions->setVisible(true);
    d->ConeSourceOptions->setEnabled(false);
    d->ConeSourceOptions->setVisible(false);
    d->SphereSourceOptions->setEnabled(false);
    d->SphereSourceOptions->setVisible(false);
  }
  else if (sourceOption == d->logic()->CONE_3D)
  {
    d->ArrowSourceOptions->setEnabled(false);
    d->ArrowSourceOptions->setVisible(false);
    d->ConeSourceOptions->setEnabled(true);
    d->ConeSourceOptions->setVisible(true);
    d->SphereSourceOptions->setEnabled(false);
    d->SphereSourceOptions->setVisible(false);  
  }
  else if (sourceOption == d->logic()->SPHERE_3D)
  {
    d->ArrowSourceOptions->setEnabled(false);
    d->ArrowSourceOptions->setVisible(false);
    d->ConeSourceOptions->setEnabled(false);
    d->ConeSourceOptions->setVisible(false);
    d->SphereSourceOptions->setEnabled(true);
    d->SphereSourceOptions->setVisible(true);
  }
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::onLogicModified()
{
  this->update();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::visualize()
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  
  // TODO: Check for empty reference volume
  if (d->InputFieldComboBox->currentNodeID() != NULL && d->OutputModelComboBox->currentNodeID() != NULL)
  {
    QProgressDialog *visualizeProgress =  new QProgressDialog(qSlicerApplication::application()->mainWindow());
    visualizeProgress->setCancelButton(0);
    visualizeProgress->setModal(true);
    visualizeProgress->setMinimumDuration(100);
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

  connect(d->InputFieldComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(inputChanged(vtkMRMLNode*)));
  connect(d->InputReferenceComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(referenceVolumeChanged(vtkMRMLNode*)));
  connect(d->OutputModelComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(outputModelChanged(vtkMRMLNode*)));

  // Glyph Parameters
  connect(d->InputGlyphPointMax, SIGNAL(valueChanged(double)), this, SLOT(setGlyphPointMax(double)));
  connect(d->InputGlyphScale, SIGNAL(valueChanged(double)), this, SLOT(setGlyphScale(double)));
  connect(d->InputGlyphThreshold, SIGNAL(valuesChanged(double, double)), this, SLOT(setGlyphThreshold(double, double)));
  connect(d->GenerateSeedButton, SIGNAL(clicked()), this, SLOT(setSeed()));
  connect(d->InputGlyphSeed, SIGNAL(valueChanged(int)), this, SLOT(setGlyphSeed(int)));  
  connect(d->GlyphSourceComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setGlyphSourceOption(int)));
  // Arrow Parameters
  connect(d->InputGlyphArrowScaleDirectional, SIGNAL(toggled(bool)), this, SLOT(setGlyphArrowScaleDirectional(bool)));
  connect(d->InputGlyphArrowScaleIsotropic, SIGNAL(toggled(bool)), this, SLOT(setGlyphArrowScaleIsotropic(bool)));  
  connect(d->InputGlyphArrowTipLength, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowTipLength(double)));
  connect(d->InputGlyphArrowTipRadius, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowTipRadius(double)));
  connect(d->InputGlyphArrowShaftRadius, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowShaftRadius(double)));  
  connect(d->InputGlyphArrowResolution, SIGNAL(valueChanged(double)), this, SLOT(setGlyphArrowResolution(double)));
  // Cone Parameters
  connect(d->InputGlyphConeScaleDirectional, SIGNAL(toggled(bool)), this, SLOT(setGlyphConeScaleDirectional(bool)));
  connect(d->InputGlyphConeScaleIsotropic, SIGNAL(toggled(bool)), this, SLOT(setGlyphConeScaleIsotropic(bool)));  
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
  //connect(d->InputContourNumber, SIGNAL(valueChanged(double)), this, SLOT(setContourNumber(double)));
  //connect(d->InputContourRange, SIGNAL(valuesChanged(double, double)), this, SLOT(setContourRange(double, double)));
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
// Set parameters
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::inputChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (pNode == NULL || this->mrmlScene() == NULL || node == NULL)
  {
    std::cerr << "Error: Unable to set input attribute" << std::endl;
    return;
  }
  
  if (node->IsA("vtkMRMLLinearTransformNode") || node->IsA("vtkMRMLBSplineTransformNode") || node->IsA("vtkMRMLGridTransformNode"))
  { 
    d->InputReferenceComboBox->setEnabled(true);
  }
  else if (node->IsA("vtkMRMLVectorVolumeNode"))
  {
    d->InputReferenceComboBox->setEnabled(false);
  }
  else
  {
    std::cerr << "Error: Unsupported input type" << std::endl;
    return;  
  }
  pNode->SetAndObserveInputNodeID(node->GetID());
  this->updateLabels();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::referenceVolumeChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene() || !node)
  {
    std::cerr << "Error: Unable to set reference volume attribute" << std::endl;
    return;
  }
  pNode->SetAndObserveReferenceVolumeNodeID(node->GetID());
  this->updateLabels();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::outputModelChanged(vtkMRMLNode* node)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene() || !node)
  {
    std::cerr << "Error: Unable to set output model attribute" << std::endl;
    return;
  }
  pNode->SetAndObserveOutputModelNodeID(node->GetID());
  this->updateLabels();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphPointMax(double pointMax)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    std::cerr << "Error: Unable to set attribute" << std::endl;
    return; 
  }
  pNode->SetGlyphPointMax(pointMax);
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
  pNode->SetGlyphSeed(rand());
  d->InputGlyphSeed->setValue(pNode->GetGlyphSeed());
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
  pNode->SetGlyphSeed(seed);
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
  pNode->SetGlyphScale(scale);
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
  pNode->SetGlyphThresholdMin(min);
  pNode->SetGlyphThresholdMax(max);
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
  pNode->SetGlyphSourceOption(option);
  this->updateGlyphSourceOptions(option);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowScaleDirectional(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphArrowScaleDirectional(state);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowScaleIsotropic(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphArrowScaleIsotropic(state);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphArrowTipLength(double length)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphArrowTipLength(length);
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
  pNode->SetGlyphArrowTipRadius(radius);
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
  pNode->SetGlyphArrowShaftRadius(radius);
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
  pNode->SetGlyphArrowResolution(resolution);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeScaleDirectional(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphConeScaleDirectional(state);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeScaleIsotropic(bool state)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphConeScaleIsotropic(state);
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModuleWidget::setGlyphConeHeight(double height)
{
  Q_D(qSlicerTransformVisualizerModuleWidget);
  vtkMRMLTransformVisualizerNode* pNode = d->logic()->GetTransformVisualizerNode();
  if (!pNode || !this->mrmlScene())
  {
    return;
  }
  pNode->SetGlyphConeHeight(height);
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
  pNode->SetGlyphConeRadius(radius);
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
  pNode->SetGlyphConeResolution(resolution);
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
  pNode->SetGlyphSphereResolution(resolution);
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
  pNode->SetGridScale(scale);
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
  pNode->SetGridSpacingMM(spacing);
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
  pNode->SetBlockScale(scale);
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
  pNode->SetBlockDisplacementCheck(state);
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
  pNode->SetContourNumber(number);
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
  //pNode->SetContourMin(min);
  //pNode->SetContourMax(max);
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
  pNode->SetContourDecimation(reduction);
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
  pNode->SetAndObserveGlyphSliceNodeID(node->GetID());
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
  pNode->SetGlyphSlicePointMax(pointMax);
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
  pNode->SetGlyphSliceThresholdMin(min);
  pNode->SetGlyphSliceThresholdMax(max);
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
  pNode->SetGlyphSliceScale(scale);
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
  pNode->SetGlyphSliceSeed(seed);
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
  pNode->SetGlyphSliceSeed(rand());
  d->InputGlyphSliceSeed->setValue(pNode->GetGlyphSliceSeed());
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
  pNode->SetAndObserveGridSliceNodeID(node->GetID());
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
  pNode->SetGridSliceScale(scale);
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
  pNode->SetGridSliceSpacingMM(spacing);
}
