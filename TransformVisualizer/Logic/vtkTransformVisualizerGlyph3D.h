#ifndef __vtkTransformVisualizerGlyph3D_h
#define __vtkTransformVisualizerGlyph3D_h

#include "vtkGlyph3D.h"
#include <vtkSmartPointer.h>
#include "vtkMinimalStandardRandomSequence.h"

#include "vtkSlicerTransformVisualizerModuleLogicExport.h"

class vtkMinimalStandardRandomSequence;

//------------------------------------------------------------------------------
class VTK_SLICER_TRANSFORMVISUALIZER_MODULE_LOGIC_EXPORT vtkTransformVisualizerGlyph3D : public vtkGlyph3D
{
public:
  vtkTypeMacro(vtkTransformVisualizerGlyph3D,vtkGlyph3D);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkTransformVisualizerGlyph3D *New();

  virtual int IsPointVisibleForListIndex(double vmag, vtkIdType ptId, int visibleListIndex);

  vtkSetMacro(PointMax,int);
  vtkGetMacro(PointMax,int);  
  vtkSetMacro(ScaleFactor,double);
  vtkGetMacro(ScaleFactor,double);
  vtkSetMacro(ScaleDirectional,bool);
  vtkGetMacro(ScaleDirectional,bool);
  vtkSetMacro(Seed,double);
  vtkGetMacro(Seed,double);
  vtkSetMacro(MagnitudeMax,double);
  vtkGetMacro(MagnitudeMax,double);
  vtkSetMacro(MagnitudeMin,double);
  vtkGetMacro(MagnitudeMin,double);
  
protected:
  vtkTransformVisualizerGlyph3D();
  ~vtkTransformVisualizerGlyph3D() {};
  
  int PointMax;
  double ScaleFactor;
  bool ScaleDirectional;
  unsigned Seed;
  double MagnitudeMax;
  double MagnitudeMin;
  
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkTransformVisualizerGlyph3D(const vtkTransformVisualizerGlyph3D&);  // Not implemented.
  void operator=(const vtkTransformVisualizerGlyph3D&);  // Not implemented.
};

#endif
