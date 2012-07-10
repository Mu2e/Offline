//
// container for the info of the extrapolated trajectory on the calorimeter
//
// $Id: TrkToCaloExtrapol.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo
//


#ifndef TrackCaloMatching_TrkToCaloExtrapol_hh
#define TrackCaloMatching_TrkToCaloExtrapol_hh

// Mu2e includes:
#include "art/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "BaBar/BbrGeom/include/BbrLorentzVectorErr.hh"

//tracker includes
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRep.hh"
#include "KalmanTrack/KalRep.hh"

//#include "BaBar/include/TrkBase/TrkDifTraj.hh"
#include "BaBar/include/TrkBase/TrkExchangePar.hh"
//#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
//#include "KalmanTests/inc/TrkRecoTrkCollection.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/TrkHelixFit.hh"

// C++ includes
#include <vector>



namespace mu2e {

typedef  art::Ptr< const TrkRecoTrk * const  >            TrkRecoTrkPtr;

struct TrkToCaloExtrapol{


private:
        int                                            _vaneId;    //vane index, runs from 0 to nVanes
        TrkRecoTrkPtr                                     _trk;
        double                             _pathLengthEntrance;
        double                                 _pathLengthExit;



public:

        TrkToCaloExtrapol():_vaneId(-1),
        _pathLengthEntrance(0.0),
        _pathLengthExit(0.0){}



        TrkToCaloExtrapol(int& vane, TrkRecoTrkPtr& trk, double& entrance, double& exit):
                _vaneId(vane),
                _trk(trk),
                _pathLengthEntrance(entrance),
                _pathLengthExit(exit){
        }
        ~TrkToCaloExtrapol(){}

        //Accessors
        int                                             vaneId() const;
        double                                            time() const;
        double                                         timeErr() const;
        double                              pathLengthEntrance() const;
        double                           pathLenghtEntranceErr() const;
        double                                  pathLengthExit() const;
        double                                              t0() const;
        double                                           t0Err() const;
        Hep3Vector                                  t0Momentum() const;
        BbrVectorErr                             t0MomentumErr() const;
        HepPoint                              entrancePosition() const;
        BbrPointErr                        entrancePositionErr() const;
        HepPoint                                  exitPosition() const;
        BbrPointErr                            exitPositionErr() const;
        Hep3Vector                                    momentum() const;
        BbrVectorErr                               momentumErr() const;
        TrkRecoTrkPtr const&                               trk() const{ return _trk; }


        void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

        bool           operator == (const TrkToCaloExtrapol & other) const ;

        bool           operator<( const TrkToCaloExtrapol& other) const{
                return ( _pathLengthEntrance< other._pathLengthEntrance);
        }

};

inline std::ostream& operator<<( std::ostream& ost,TrkToCaloExtrapol const& t){
        t.print(ost,false);
        return ost;
}

} // namespace mu2e

#endif /* TrackCaloMatching_TrkToCaloExtrapol_hh */
