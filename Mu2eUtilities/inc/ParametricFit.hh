#ifndef _MU2E_UTILITIES_PARAMETRICFIT_HH
#define _MU2E_UTILITIES_PARAMETRICFIT_HH

#include <vector>

#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "Math/VectorUtil.h"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "TMatrix.h"

namespace mu2e { struct ComboHit; }

using namespace mu2e;
namespace ParametricFit{
       
	XYZVec pointOnLineFromX(XYZVec lineStartPoint, XYZVec lineEndPoint, double x,XYZVec outputPoint);
	bool LineToLineCA(XYZVec& FirstLinePos, XYZVec& FirstLineDir, 
	 XYZVec& SecondLinePos, XYZVec& SecondLineDir, 
	 XYZVec& closestPointOnFirstLine, XYZVec& closestPointOnSecondLine);
	double LineToLineDCA(XYZVec& firstLineStartPoint, XYZVec& firstLineEndPoint,XYZVec& secondLineStartPoint, XYZVec& secondLineEndPoint, double& dca);
	
        
	std::vector<XYZVec>GetAxes(XYZVec TrackDirection);
	TrackAxes GetTrackAxes(XYZVec TrackDirection);
	XYZVec GetXPrime(XYZVec track_dir);
	XYZVec GetYPrime(XYZVec OrthX, XYZVec YPrime);
	XYZVec GetXDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	XYZVec GetYDoublePrime(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	
	void TestConditions(XYZVec XPrime, XYZVec YPrime, XYZVec ZPrime);
	XYZVec ConvertToPrimePoint(ComboHit* chit, TrackAxes axes);
	
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
        
	}

#endif
