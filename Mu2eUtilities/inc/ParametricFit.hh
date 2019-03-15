#ifndef _MU2E_UTILITIES_PARAMETRICFIT_HH
#define _MU2E_UTILITIES_PARAMETRICFIT_HH
// Fitting Routines for parametric line
// Author: S. Middleton
// Date: March 2019

//ROOT

#include "Math/VectorUtil.h"
#include "TMatrix.h"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/DriftCircle.hh"
using namespace mu2e;
namespace ParametricFit{

	double GettMin(XYZVec& point, XYZVec& starting_point, XYZVec& end_point);
	XYZVec PointToLineCA(XYZVec& point, XYZVec& starting_point, XYZVec& end_point);
	double PointToLineDCA(XYZVec& point, XYZVec& line_starting_point, XYZVec& line_end_point);
	XYZVec ParellelVector(XYZVec lineStartPoint, XYZVec lineEndPoint);

	//CosmicTrack ConstructTrack();
	XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint);

	bool LineToLineCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint, 
	  XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, 
	  XYZVec& closestPointOnFirstLine, XYZVec& closestPointOnSecondLine);

	double LineToLineDCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint,XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, double& dca);
        XYZVec MajorAxis(ComboHit* Hit, XYZVec track_dir);

	XYZVec MinorAxis(ComboHit* Hit, XYZVec track_dir);

	double HitErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir);

        double HitErrorY(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir);
        
        double TotalHitError(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec track_dir);
        
        int GetDOCASign(XYZVec track_dir, XYZVec point);
       
        double GetResidualX(double A0, double A1, XYZVec track_dir, XYZVec point);
        double GetResidualY( double B0, double B1, XYZVec track_dir, XYZVec point);
        
	double GetResidualError(XYZVec Major_Axis, XYZVec Minor_Axis, XYZVec track_direction, XYZVec& POCA, double DOCA);

	double Get2DParameter(int i, TMatrixD alpha);	

	TMatrixD GetGamma(double G00, double Gij, double G11);

	TMatrixD GetBeta(double	b0,double b1);

	}

#endif
