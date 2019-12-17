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
struct FitResult{
        public:
		std::vector<std::string> names;
		std::vector<double> bestfit;
		std::vector<double> bestfiterrors;
		std::vector<double> bestfitcov;

		std::vector<double> StartDOCAs;
		std::vector<double> StartTimeResiduals;

		std::vector<double> GaussianEndDOCAs;
		std::vector<double> GaussianEndTimeResiduals;

		std::vector<double> FullFitEndDOCAs;
		std::vector<double> FullFitEndTimeResiduals;

		std::vector<double> RecoAmbigs;
		
		double NLL;
	
};

namespace MinuitDriftFitter {
	
	FitResult DoFit(int diag, CosmicTrackSeed trackseed, StrawResponse srep, double doca_cut, unsigned int MinNCh_cut, int LogLcut, double _gaussTres, double _maxTres);

}


#endif
