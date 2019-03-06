#ifndef _MU2E_UTILITIES_PARAMETRICFIT_HH
#define _MU2E_UTILITIES_PARAMETRICFIT_HH
// Fitting Routines for parametric line
// Author: S. Middleton
// Date: March 2019

//ROOT

#include "Math/VectorUtil.h"
#include "RecoDataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/DriftCircle.hh"
using namespace mu2e;
namespace ParametricFit{

	double GettMin(XYZVec& point, XYZVec& starting_point, XYZVec& end_point);
	void PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point, XYZVec& closestPointOnLine, bool finite);
	double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point, double& dca, bool finiteLine);
	XYZVec ParellelVector(XYZVec lineStartPoint, XYZVec lineEndPoint);

	//CosmicTrack ConstructTrack();
	XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint,bool finiteLine);

	bool LineToLineCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint, 
	  XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, 
	  XYZVec& closestPointOnFirstLine, XYZVec& closestPointOnSecondLine, bool finiteLine);

	double LineToLineDCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint,XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, double& dca, bool finiteLine);
       
	vector<mu2e::StraightTrack*> Calculate2DLineFits(std::vector<mu2e::DriftCircle> circles, double pValCut, long LR);
        void Construct3DTrack(StraightTrack* xyLineFit, StraightTrack* zrLineFit, CosmicTrack* track);

	}

#endif
