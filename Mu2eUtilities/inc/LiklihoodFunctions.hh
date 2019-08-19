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
		std::vector<double> StartDOCAs;
		std::vector<double> StartTimeResiduals;
		std::vector<double> EndDOCAs;
		std::vector<double> EndTimeResiduals;
		
	
};

namespace LiklihoodFunctions {
	
	EndResult DoFit(CosmicTrackSeed trackseed, StrawResponse srep);

}


#endif
