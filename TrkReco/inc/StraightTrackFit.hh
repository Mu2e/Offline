//
// Class to perform Straight fit
// Original author: Sophie Middleton 
//
// $Id: StraightTrackFit.hh 
// $Author: sophie $ 
// $Date: 2018/01/12 18:56:10 $
//
#ifndef TrkReco_StraightTrackFit_HH
#define TrkReco_StraightTrackFit_HH

// framework
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"

// Mu2e objects
#include "TrkReco/inc/StraightTrackFinderData.hh"
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
namespace mu2e 
{
    class Tracker;
    class StraightTrackFit
    {
     public:
                explicit StraightTrackFit(fhicl::ParameterSet const&);
    		virtual ~StraightTrackFit();

                bool initStraightTrack(StraightTrackFinderData& TrackData);
                StraightTrack* InitLine(const ComboHit *FirstP1, const ComboHit *LastP2,StraightTrack* line);
                //Step 1: Begin Fit- initializes the fit routine:
                void BeginFit(StraightTrackFinderData& TrackData);
                
                //Step 2: RunFitChi2-holds the functions to find initial line, update, refine and add in drift
                void RunFitChi2(StraightTrackFinderData& trackData);
		//Step 3: Fit All - finds the chi-squared anf line information when all hits in track taken into account. This will be the initial chi-squared value.
		StraightTrack* FitAll(StraightTrackFinderData& trackData,StraightTrack* track, int WeightMode);
		
		void MulitpleTrackResolver(StraightTrackFinderData& trackData,StraightTrack* track);

		std::vector<double> Optimize(StraightTrackFinderData& trackData, StraightTrack* track, std::vector<std::vector<double>> all_hit_list);
                
                void UpdateFitErrors(std::vector<double> x, std::vector<double> y,std::vector<double> z, std::vector<double> err, StraightTrack* track,TMatrixD cov_x, std::vector<XYZVec> maj,std::vector<XYZVec> min);

                bool goodTrack(StraightTrack* track);
                float pointToLineDCA(StraightTrack* track, StrawHit hit);
		void DriftCorrection(StraightTrackFinderData& trackData);
		StraightTrack* Hough2D(StraightTrackFinderData& trackData);

                float evalWeightXY  (const ComboHit& Hit, const StraightTrack& track);
                const Tracker*            _tracker;
    		void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
                
		
	private:
		
                //list function:
               
  		bool use(ComboHit const&) const;
    		void setOutlier(ComboHit&) const;
                float hitWeight(ComboHit const& hhit) const;
                int _dim; //dimensions of fit
		int _diag;
    		int _debug;		    // debug level
                StrawHitFlag _useflag, _dontuseflag;
	        int _minresid;
    		int      _minnsh;  // minimum # of StrawHits
    		unsigned _minCHHits; // minimum # of CH hits - should be at least 2 for fit to work....
		unsigned _n_outliers; //number of significant outliers/number of hits in track....helps with multiplicity(?)
    		float _maxresid; // max allowed pull which a hit can have to be classed as "good"
    		unsigned _maxniter; // maxium # of iterations to global minimum   
                
    		float _maxdxy; // maximum distance in hits after fit
    		float _maxchi2; //maximum allowed chi2
		
    
  };//end Fit class

}
#endif
