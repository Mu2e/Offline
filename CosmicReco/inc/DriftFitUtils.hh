#ifndef _COSMIC_RECO_DriftFitUtils_HH
#define _COSMIC_RECO_DriftFitUtils_HH
//Date: August 2019
//Author: S. Middletin
//Purpose: Add in drift fit functionality in seperate file for ease
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/TrackerConditions/inc/StrawDrift.hh"

//Tracker Drift Conditions:
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/StrawPhysics.hh"
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
//Utilities:
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"


using namespace mu2e;

namespace DriftFitUtils{
          TrackEquation ConvertFitToDetectorFrame(TrackAxes axes, XYZVectorF Position, XYZVectorF Direction);
        double GetTestDOCA(ComboHit const& chit,XYZVectorF const& pos, XYZVectorF const& dir,  const Tracker* tracker);
        double GetTestDOCA(ComboHit const& chit,double a0, double a1, double b0, double b1, const Tracker* tracker);
        double GetRPerp(StrawResponse const& srep, ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker);
        double GetDriftDistance(StrawResponse const& srep, ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker);
        int GetAmbig(ComboHit const& chit, XYZVectorF const& pos, XYZVectorF const& dir,  const Tracker* tracker);
        int GetAmbig(ComboHit const& chit, double a0, double a1, double b0, double b1,  const Tracker* tracker);
          double GetPropVelocity(StrawResponse const& srep, ComboHit const& chit);
        double GetPropTime(ComboHit const& chit, double vprop, const Tracker* tracker);
          double TimeResidualTrans(double doca);
          double TimeResidualLong(double doca, StrawResponse const& srep, double t0, ComboHit const& chit,  const Tracker* tracker);
          double TimeResidual(double doca, StrawResponse const& srep, double t0, ComboHit const& hit,  const Tracker* tracker);


 }

#endif
