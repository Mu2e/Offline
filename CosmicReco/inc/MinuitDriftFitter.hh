#ifndef _COSMIC_RECO_MINUITDRIFTFITTER_HH
#define _COSMIC_RECO_MINUITDDRIFTFITTER_HH

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
	FitResult DoFit(int diag, CosmicTrackSeed& tseed, StrawResponse const& srep, const Tracker* tracker, double doca_cut, unsigned int MinNCh_cut, int LogLcut, double _gaussTres, double _maxTres);
        void DoDriftTimeFit(int diag, CosmicTrackSeed& tseed, StrawResponse const& srep, const Tracker* tracker );

    }


#endif


