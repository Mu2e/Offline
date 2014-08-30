#ifndef RecoDataProducts_TrackSeed_hh
#define RecoDataProducts_TrackSeed_hh
//
// out data of the time pattern recognition algorithm to seed the Kalman filter (used for the ITracker)
//
// $Id: TrackSeed.hh,v 1.2 2014/08/30 12:16:39 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:16:39 $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>

// Mu2e includes
#include "art/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/HelixVal.hh"

namespace mu2e {

typedef art::Ptr<StrawHit> StrawHitPtr;
typedef art::Ptr<TrackerHitTimeCluster> TrackerHitTimeClusterPtr;

struct TrackSeed{

        std::vector<StrawHitPtr> _selectedTrackerHits;
        TrackerHitTimeClusterPtr _relatedTimeCluster;

        HelixVal _fullTrkSeed;
        std::vector<mu2e::HelixVal> _loopSeeds;
        double _t0;
        double _errt0;


public:

        TrackSeed():
                _t0(0.),
                _errt0(0.) {
        }
        double d0() const {return _fullTrkSeed._d0;}
        double phi0() const {return _fullTrkSeed._phi0;}
        double omega() const  {return _fullTrkSeed._omega;}
        double z0() const  {return _fullTrkSeed._z0;}
        double tanDip() const  {return _fullTrkSeed._tanDip;}

        double t0() const  {return _t0;}
        double errt0() const  {return _errt0;}
        size_t nLoops() const {return _loopSeeds.size(); }

        // Print contents of the object.
        void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
};

inline std::ostream& operator<<( std::ostream& ost,
                TrackSeed const& seed){
        ost<<"Reconstructed track with parameters: d0= " << seed.d0() << " phi0= "<< seed.phi0() << " omega= "<< seed.omega() << " z0 = "<< seed.z0() << " tanDip= "<< seed.tanDip() << std::endl;
        ost<<"Cov Matrix:"<<std::endl;
        for (int row=0; row<5;++row) {
                for (int col=0;col<5;++col) {
                        ost<<seed._fullTrkSeed._covMtrx[row][col]<<" ";
                }
                ost<<std::endl;
        }
        ost<<"with "<<seed.nLoops()<<" loops and potential t0 "<<seed.t0()<<" with error "<<seed.errt0()<<std::endl;
        return ost;
}


} // namespace mu2e

#endif /* RecoDataProducts_TrackSeed_hh */
