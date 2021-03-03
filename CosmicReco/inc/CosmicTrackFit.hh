//Author: S Middleton
//Purpose: Fit cosmic tracks within the tracker

#ifndef CosmicReco_CosmicTrackFit_HH
#define CosmicReco_CosmicTrackFit_HH

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

//Mu2e Cosmics:
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

// Products
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"

//Drift:
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawDrift.hh"

// Math
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"

//C++
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>

// Framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

//ROOT
#include "TMatrixD.h"

namespace mu2e 
{
    class Tracker;
    class CosmicTrackFit
    {
     public:	

	struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> Npara{Name("NParameters"),Comment("number of fit parameters used")};
	      fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),0};
	      fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints"),0};
	      fhicl::Atom<string> dontuseflag {Name("DoNotUseFlag"),Comment("if set to OK then save the track")};
              fhicl::Atom<unsigned> minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits ")};
    	      fhicl::Atom<unsigned> minnch {Name("minNComboHits"), Comment("number of combohits allowed")};
	      fhicl::Atom<unsigned> n_outliers{Name("Noutliers"),Comment("maximum number of outliers allowed in track fit")};
    	      fhicl::Atom<unsigned> maxniter{Name("maxNiter"), Comment("Maximum allowed number of iterations before considered unconverged in seed fit")};
    	      
              fhicl::Atom<unsigned> minNHitsTimeCluster{Name("minNHitsTimeCluster"),Comment("minium allowed time cluster")};
    	      fhicl::Atom<float> max_seed_chi2{Name("MaxSeedChi2DOF"),Comment("maximum chi 2/dof for seed")};
	      fhicl::Atom<float> max_chi2_change{Name("MaxDeltaChi2"),Comment("The maxiumum allowed change in chi2 before convergeing seed fit")};
	      fhicl::Atom<float> max_position_deviation{Name("MaxPosDev"),Comment("The maxiumum allowed change in position correlated parameters between seed fit iterations ") };
	      fhicl::Atom<float> maxHitDOCA{Name("MaxDOCA"),Comment("The maxiumum allowed DOCA to wire for any hit used for full drift fit") };
	      fhicl::Atom<float> maxLogL{Name("MaxLogL"),Comment("The maxiumum allowed outcome of Minuit fit routine") };
	      fhicl::Atom<float> gaussTres{Name("GaussianSeedTimeResolution"),Comment("The resolution of the Gaussian seed fit in time") };
              fhicl::Atom<float> maxTres{Name("MaxTimeResidual"),Comment("The maxiumum allowed time residual for any hit used for full drift fit") };
	      fhicl::Atom<float> maxd{Name("MaxTrackLength"),Comment("The maxiumum allowed length of track")};
	      fhicl::Atom<float> maxpull{Name("MaxHitPullForSeed"),Comment("The maxiumum allowed combo hit pull from fit")};
              fhicl::Atom<bool> UseTSeedDirection{Name("UseTSeedDirection"),Comment("Uses the direction in the input tseed to initialize fit")};
    	};
		
		explicit CosmicTrackFit(const Config& conf);
    		virtual ~CosmicTrackFit(){};

                bool initCosmicTrack(const char* title, CosmicTrackSeed& tseed, ComboHitCollection &combohits);
		std::vector<XYZVec> SortPoints(std::vector<XYZVec> pointY);
                XYZVec InitLineDirection(const ComboHit *ch0, const ComboHit *chN);
                
                XYZVec LineDirection(double a0, double a1, const ComboHit *ch0, const ComboHit *chN, XYZVec ZPrime);
                XYZVec ConvertPointToDetFrame(XYZVec vec);

                XYZVec GetTrackDirection(std::vector<XYZVec> hitXYZ, XYZVec XDoublePrime, XYZVec YDoublePrime, XYZVec ZPrime); 
                void BeginFit(const char* title, CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs);
                void RunFitChi2(const char* title, CosmicTrackSeed& tseed, ComboHitCollection &combohits);
                void FitAll(const char* title, CosmicTrackSeed &tseed, ComboHitCollection &combohits, CosmicTrack* cosmictrack);
                void FillTrackHitCollections(CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs);

		void ConvertFitToDetectorFrame(TrackAxes axes, XYZVec Position, XYZVec Direction, CosmicTrack* cosmictrack, bool isseed, bool det);
		
                bool goodTrack(CosmicTrack& track);
		void DriftFit(CosmicTrackSeed& tseed, StrawResponse const& srep);
		
                const Tracker*            _tracker;
    		void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
                
		bool use_hit(ComboHit const&) const;
  		bool use_track(double length) const;
	private:
		Config _conf;
  		
                unsigned _Npara;
		int _diag;
    		int _debug;		  
                StrawHitFlag _useflag, _dontuseflag; //some flags for removn
    		unsigned _minnsh;  // minimum # of StrawHits
    		unsigned _minnch; // minimum # of CH hits - should be at least 2 for fit to work....
		unsigned _n_outliers; //number of significant outliers/number of hits in track..
    		unsigned _maxniter; // maxium # of iterations to global minimum         
 	        float _max_seed_chi2;// max chi2 total of seed track
		float _max_chi2_change; // once we are lower than this we can say its converged
	        float _max_position_deviation;//maximum hcange from initial choice of a0 and b0 (in mm)
	        float _maxHitDOCA ;//Max allowed DOCA for gaussian seed
		float _maxLogL; //Maximum allowed liklihood in drift model
		float _gaussTres;//resolution for intial seed
		float _maxTres;//maximum allowed time residual in drift fit
		float _maxd;//unused
		float _maxpull;
                bool _useTSeedDirection;
  };
	

}
#endif
