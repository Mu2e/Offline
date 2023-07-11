#ifndef _MU2E_UTILITIES_PARAMETRICFIT_HH
#define _MU2E_UTILITIES_PARAMETRICFIT_HH

#include <vector>

#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Math/VectorUtil.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "TMatrix.h"

namespace mu2e { struct ComboHit; }

using namespace mu2e;
namespace ParametricFit{

        XYZVectorF pointOnLineFromX(XYZVectorF lineStartPoint, XYZVectorF lineEndPoint, double x,XYZVectorF outputPoint);
        bool LineToLineCA(XYZVectorF& FirstLinePos, XYZVectorF& FirstLineDir,
         XYZVectorF& SecondLinePos, XYZVectorF& SecondLineDir,
         XYZVectorF& closestPointOnFirstLine, XYZVectorF& closestPointOnSecondLine);
        double LineToLineDCA(XYZVectorF& firstLineStartPoint, XYZVectorF& firstLineEndPoint,XYZVectorF& secondLineStartPoint, XYZVectorF& secondLineEndPoint, double& dca);


        std::vector<XYZVectorF>GetAxes(XYZVectorF TrackDirection);
        TrackAxes GetTrackAxes(XYZVectorF TrackDirection);
        XYZVectorF GetXPrime(XYZVectorF track_dir);
        XYZVectorF GetYPrime(XYZVectorF OrthX, XYZVectorF YPrime);
        XYZVectorF GetXDoublePrime(XYZVectorF XPrime, XYZVectorF YPrime, XYZVectorF ZPrime);
        XYZVectorF GetYDoublePrime(XYZVectorF XPrime, XYZVectorF YPrime, XYZVectorF ZPrime);

        void TestConditions(XYZVectorF XPrime, XYZVectorF YPrime, XYZVectorF ZPrime);
        XYZVectorF ConvertToPrimePoint(ComboHit* chit, TrackAxes axes);

        std::vector<double> GetErrors(ComboHit* Hit, XYZVectorF XAxis, XYZVectorF YAxis);
        XYZVectorF MajorAxis(ComboHit* Hit);
        XYZVectorF MinorAxis(ComboHit* Hit);
        double HitErrorX(ComboHit* Hit, XYZVectorF major_axis, XYZVectorF minor_axis, XYZVectorF XPrime);
        double HitErrorY(ComboHit* Hit, XYZVectorF major_axis, XYZVectorF minor_axis, XYZVectorF YPrime);
        double TotalHitError(ComboHit* Hit, XYZVectorF major_axis, XYZVectorF minor_axis, XYZVectorF XPrime, XYZVectorF YPrime);

        double GetHitChi2(double A0, double A1, double errorX, XYZVectorF point, XYZVectorF XDoublePrime, XYZVectorF ZPrime);
        double GetGlobalChi2(double a0, double a1, double b0, double b1, XYZVectorF prime_point, double errX, double errY, int DOF);
        double GetResidualX(double A0, double A1,  XYZVectorF point);
        double GetResidualY( double B0, double B1,  XYZVectorF point);
        double GetTotalResidual(double resid_x,double resid_y);
        double GetResidualError(XYZVectorF Major_Axis, XYZVectorF Minor_Axis, XYZVectorF track_direction, XYZVectorF& POCA, double DOCA);

        }

#endif
