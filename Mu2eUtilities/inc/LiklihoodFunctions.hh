#ifndef _MU2E_UTILITIES_LIKLIHOODFUNCTIONS_HH
#define _MU2E_UTILITIES_LIKLIHOODFUNCTIONS_HH
// Author: S. Middleton 
// Date: July 2019
//Purpose: Will pass PDF function to Minuit 
#include "TrackerConditions/inc/StrawDrift.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
//For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
//ROOT
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
//Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;
struct EndResult{
        public:
		std::vector<std::string> names;
		std::vector<double> bestfit;
		std::vector<double> bestfiterrors;
	
};

namespace LiklihoodFunctions {
	//double TimeResiduals(double drifttime, double doca, double time, double time_offset);
        //double GetDOCA(TrkPoca poca);
       // TrkPoca GetPOCA(Straw const&  straw, std::vector<XYZVec> TrackAxes, XYZVec track_position, XYZVec track_direction);
        //double GetPhi(Straw const&  straw, XYZVec track_direction, std::vector<XYZVec> TrackAxes, double doca);
	EndResult DoFit(CosmicTrackSeed trackseed, std::vector<double> times,  std::vector<XYZVec> position, std::vector<double> errorsX, std::vector<double> errorsY, std::vector<Straw> straws);

}


#endif
