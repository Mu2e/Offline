#ifndef _MU2E_UTILITIES_DriftFitUtils_HH
#define _MU2E_UTILITIES_DriftFitUtils_HH
//Date: August 2019
//Author: S. Middletin
//Purpose: Add in drift fit functionality in seperate file for ease
#include "TrackerConditions/inc/StrawDrift.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
//Tracker Drift Conditions:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
//Utilities:
#include "Mu2eUtilities/inc/ParametricFit.hh"
//For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"


using namespace mu2e;

namespace DriftFitUtils{
  	TrackEquation ConvertFitToDetectorFrame(TrackAxes axes, XYZVec Position, XYZVec Direction);
	double GetTestDOCA(ComboHit chit,XYZVec pos, XYZVec dir);
	double GetTestDOCA(ComboHit chit,double a0, double a1, double b0, double b1);
        int GetAmbig(ComboHit chit, XYZVec pos, XYZVec dir);
        int GetAmbig(ComboHit chit, double a0, double a1, double b0, double b1);
  	double GetPropVelocity(StrawResponse::cptr_t srep, ComboHit chit); 
	double GetPropTime(ComboHit chit, double vprop);
  	double TimeResidualTrans(double doca);
  	double TimeResidualLong(double doca, StrawResponse::cptr_t srep, double t0, ComboHit chit);
  	double TimeResidual(double doca, StrawResponse::cptr_t srep, double t0, ComboHit hit);
  	
  
 }

#endif
