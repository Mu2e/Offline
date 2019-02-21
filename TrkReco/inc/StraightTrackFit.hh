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
    class TTracker;
    class StraightTrackFit
    {
     public:
                explicit StraightTrackFit(fhicl::ParameterSet const&);
    		virtual ~StraightTrackFit();

                bool initStraightTrack(StraightTrackFinderData& TrackData);
                StraightTrack* InitLine(const ComboHit *FirstP1, const ComboHit *LastP2,StraightTrack* line);
                void BeginFit(StraightTrackFinderData& TrackData);
                
                void Fit(StraightTrackFinderData& trackData, double chi2); 
                void FitChi2(StraightTrackFinderData& trackData);
                StraightTrack* refineFit(StraightTrackFinderData& trackData,StraightTrack* track, int WeightMode);
                
                void UpdateFitErrors(std::vector<double> x, std::vector<double> y,std::vector<double> z, std::vector<double> err, StraightTrack* track,TMatrixD cov_x, std::vector<XYZVec> maj,std::vector<XYZVec> min);
                bool goodTrack(StraightTrack* track);
               
		void add_drift(StraightTrackFinderData& trackData);
                float evalWeightXY  (const ComboHit& Hit, const StraightTrack& track);
                const TTracker*            _tracker;
    		void  setTracker    (const TTracker*    Tracker) { _tracker     = Tracker; }
                
		
	private:
		
                //list function:
               
  		bool use(ComboHit const&) const;
    		void setOutlier(ComboHit&) const;
                float hitWeight(ComboHit const& hhit) const;
                int _dim; //dimensions of fit
		int _diag;
    		int _debug;		    // debug level
                StrawHitFlag _useflag, _dontuseflag;
    		int      _minnsh;  // minimum # of StrawHits
    		//unsigned _minnhit; // minimum # of hits to work with
    		float _minxyresid; // minimum distance used in the track fit to be clusterized. units are mm
    		unsigned _maxniter; // maxium # of iterations to global minimum   
                //float _minzsep, _maxzsep; // Z separation of points for pitch estimate
    		float _maxdxy; // maximum distance in hits after fit
    		float _maxchi2xy;
		
    
  };//end Fit class

}
#endif
