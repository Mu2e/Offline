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
	int GetDOCASign(XYZVec track_dir, XYZVec point);

	std::vector<XYZVec>GetAxes(XYZVec TrackDirection);
	XYZVec GetXPrime(XYZVec track_dir);
	XYZVec GetYPrime(XYZVec OrthX, XYZVec YPrime);
	XYZVec GetXDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	XYZVec GetYDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	void TestConditions(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	
	std::vector<double> GetErrors(ComboHit* Hit, XYZVec XAxis, XYZVec YAxis);
	XYZVec MajorAxis(ComboHit* Hit);
	XYZVec MinorAxis(ComboHit* Hit);
	double HitErrorX(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec XPrime);
        double HitErrorY(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec YPrime);
        double TotalHitError(ComboHit* Hit, XYZVec major_axis, XYZVec minor_axis, XYZVec XPrime, XYZVec YPrime);
        
        double GetHitChi2(double A0, double A1, double errorX, XYZVec point, XYZVec XDoublePrime, XYZVec ZPrime);
        double GetGlobalChi2(double a0, double a1, double b0, double b1, XYZVec prime_point, double errX, double errY, int DOF);
        double GetResidualX(double A0, double A1,  XYZVec point);
        double GetResidualY( double B0, double B1,  XYZVec point);
        double GetTotalResidual(double resid_x,double resid_y);
	double GetResidualError(XYZVec Major_Axis, XYZVec Minor_Axis, XYZVec track_direction, XYZVec& POCA, double DOCA);
        
        double GetMCResidualX(double A0, double A1, XYZVec MCTrackDirection, XYZVec X, XYZVec point);
        double GetMCResidualY(double A0, double A1, XYZVec MCTrackDirection, XYZVec Y, XYZVec point);
	

	}

#endif
