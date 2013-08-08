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
#include <QtPlugin>

// TransformVisualizer Logic includes
#include <vtkSlicerTransformVisualizerLogic.h>

// TransformVisualizer includes
#include "qSlicerTransformVisualizerModule.h"
#include "qSlicerTransformVisualizerModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerTransformVisualizerModule, qSlicerTransformVisualizerModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_TransformVisualizer
class qSlicerTransformVisualizerModulePrivate
{
public:
  qSlicerTransformVisualizerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerTransformVisualizerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModulePrivate
::qSlicerTransformVisualizerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerTransformVisualizerModule methods

//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModule
::qSlicerTransformVisualizerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerTransformVisualizerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerTransformVisualizerModule::~qSlicerTransformVisualizerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerTransformVisualizerModule::helpText()const
{
  return "The Transform Visualizer module visualizes transforms using various options. The module can visualize any transform (linear transform, B-spline deformable transform, any other non-linear transform) or vector volume. See <a>http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/TransformVisualizer</a> for more information.";
}

//-----------------------------------------------------------------------------
QString qSlicerTransformVisualizerModule::acknowledgementText()const
{
  return "This work is part of SparKit project, funded by Cancer Care Ontario (CCO)'s ACRU program and Ontario Consortium for Adaptive Interventions in Radiation Oncology (OCAIRO).";
}

//-----------------------------------------------------------------------------
QStringList qSlicerTransformVisualizerModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Franklin King") << QString("Andras Lasso") << QString("Csaba Pinter");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerTransformVisualizerModule::icon()const
{
  return QIcon(":/Icons/TransformVisualizer.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerTransformVisualizerModule::categories() const
{
  return QStringList() << "Registration";
}

//-----------------------------------------------------------------------------
QStringList qSlicerTransformVisualizerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerTransformVisualizerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerTransformVisualizerModule
::createWidgetRepresentation()
{
  return new qSlicerTransformVisualizerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerTransformVisualizerModule::createLogic()
{
  return vtkSlicerTransformVisualizerLogic::New();
}
