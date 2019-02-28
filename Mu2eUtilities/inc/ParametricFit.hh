#ifndef _MU2E_UTILITIES_PARAMETRICFIT_HH
#define _MU2E_UTILITIES_PARAMETRICFIT_HH
// Fitting Routines for parametric line
// Author: S. Middleton
// Date: March 2019

//ROOT

#include "Math/VectorUtil.h"
#include "RecoDataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
//using namespace mu2e;

namespace ParametricFit{
 
void PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point, XYZVec& closestPointOnLine, bool finite);
double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point, double& dca, bool finiteLine);
XYZVec ParellelVector(XYZVec lineStartPoint, XYZVec lineEndPoint);
double GetResidual();
//CosmicTrack ConstructTrack();
XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint,bool finiteLine);


}

#endif
