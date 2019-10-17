
#ifndef TrkReco_CosmicTrackFit_HH
#define TrkReco_CosmicTrackFit_HH

// framework
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

//Mu2e Cosmics:
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"

// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
//Drift:
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
// Mu2e objects
#include "TrkReco/inc/CosmicTrackFinderData.hh"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"

//C++
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>
#include <algorithm>

//ROOT
#include "TMatrixD.h"
#include "TH2F.h"
#include "TCanvas.h"

namespace mu2e 
{
    class Tracker;
    class CosmicTrackFit
    {
     public:
                explicit CosmicTrackFit(fhicl::ParameterSet const&);
    		virtual ~CosmicTrackFit();
                bool initCosmicTrack(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics);
                XYZVec InitLineDirection(const ComboHit *ch0, const ComboHit *chN);
                XYZVec InitLineDirection( StrawDigiMC const& ch0,  StrawDigiMC const& chN, XYZVec reco, bool is_prime) ;
                XYZVec LineDirection(double a0, double a1, const ComboHit *ch0, const ComboHit *chN, XYZVec ZPrime);
                XYZVec GetTrackDirection(std::vector<XYZVec> hitXYZ, XYZVec XDoublePrime, XYZVec YDoublePrime, XYZVec ZPrime); 
                //Step 1: Begin Fit- initializes the fit routine:
                void BeginFit(const char* title, CosmicTrackFinderData& TrackData, CosmicTrackFinderTypes::Data_t& diagnostics);
                //Step 2: RunFitChi2-holds the functions to find initial line, update, refine 
                void RunFitChi2(const char* title, CosmicTrackFinderData& trackData, CosmicTrackFinderTypes::Data_t& diagnostics);
                std::vector<XYZVec> SortPoints(std::vector<XYZVec> pointY);
		//Step 3: Fit All - finds the chi-squared anf line information when all hits in track taken into account. This will be the initial chi-squared value.
		void FitAll(const char* title, CosmicTrackFinderData& trackData,CosmicTrack* track, CosmicTrackFinderTypes::Data_t& diagnostics);
		//Step 4: Do the Chi2 fitting
		void FitXYZ(CosmicTrackFinderData& trackData,CosmicTrack* track, CosmicTrackFinderTypes::Data_t& diagnostics);
		void ConvertFitToDetectorFrame(CosmicTrackFinderData& trackData, TrackAxes axes, XYZVec Position, XYZVec Direction, CosmicTrack* cosmictrack, bool isseed, bool det);
		//Step 5: validation of algorithm -  some functions to help
		float PDF(float chisq, float ndf);
		float chi_sum(float chisq, float ndf);
		float CDF(float chisq, float ndf);
		
		//Some functions to extract MC truth
		void FitMC(CosmicTrackFinderData& trackData, CosmicTrack* cosmictrack, bool Det);
		void TransformMC(CosmicTrackFinderData& trackData, TrackAxes Axes, CosmicTrack* cosmictrack, bool is_seed);
                bool goodTrack(CosmicTrack* track);
                
                //Step 6: Drift fit - calls to logL minuit code utility
		void DriftFit(CosmicTrackFinderData& trackData);
		
                
                const Tracker*            _tracker;
    		void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
                StrawResponse _srep;
                void  setStrawResponse (StrawResponse rep) {_srep = rep;}
	private:
		
  		bool use_hit(ComboHit const&) const;
  		bool use_track(double length) const;
    		void setOutlier(ComboHit&) const;
                float hitWeight(ComboHit const& hhit) const;
                unsigned _Npara;
		int _diag;
		int _mcdiag;
    		int _debug;		    // debug level
                StrawHitFlag _useflag, _dontuseflag;
    		unsigned _minnsh;  // minimum # of StrawHits
    		unsigned _minCHHits; // minimum # of CH hits - should be at least 2 for fit to work....
		unsigned _n_outliers; //number of significant outliers/number of hits in track....helps with multiplicity(?)
    		unsigned _maxniter; // maxium # of iterations to global minimum   
		float _maxpull; //maximum allowed hit pull (res/reserror)             
    		float _maxd; // maximum distance in hits to begin fit
		
    		float _maxchi2; //maximum allowed chi2
		float _max_chi2_change; // once we are lower than this we can say its converged
	        float _max_seed_chi2;// max chi2 total of seed track
    		float _max_position_deviation;//maximum hcange from initial choice of a0 and b0 (in mm)
	        float _maxHitDOCA ;//Max allowed DOCA for gaussian seed
		float _minMomTrue;//Minimum true momentum (unused)
		int _maxLogL; //Maximum allowed liklihood in drift model
		unsigned _minCHStrawFull; //Minimum number of hits remaining after Gaussian seed
		double _gaussTres;//resolution for intial seed
		double _maxTres;//maximum allowed time residual in drift fit
  };//end Fit class
	

}
#endif
