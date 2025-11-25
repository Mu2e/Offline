//
// container for the info of the extrapolated trajectory on the calorimeter
//
//
// Original author G. Pezzullo
//


#ifndef TrackCaloMatching_TrkToCaloExtrapol_hh
#define TrackCaloMatching_TrkToCaloExtrapol_hh

// Mu2e includes:
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/BbrPointErr.hh"
#include "BTrk/BbrGeom/BbrLorentzVectorErr.hh"

//tracker includes
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

// C++ includes
#include <vector>

namespace mu2e {

  //  typedef  art::Ptr<const KalRep * const>  EKalRepPtr;
  //  typedef  art::Ptr<KalRep>  EKalRepPtr;

  struct TrkToCaloExtrapol{

  private:
    int                                    _diskId;
    int                                    _trackNumber;       // track numeber
    KalRepPtr                              _trk;
    double                                 _pathLengthEntrance;
    double                                 _pathLengthExit;



  public:

    TrkToCaloExtrapol():_diskId(-1),
                        _pathLengthEntrance(0.0),
                        _pathLengthExit(0.0){}



    TrkToCaloExtrapol(int disk, int trkNumber,
                      KalRepPtr& trk, double entrance, double exit):
      _diskId(disk),
      _trackNumber(trkNumber),
      _trk(trk),
      _pathLengthEntrance(entrance),
      _pathLengthExit(exit){
    }
    ~TrkToCaloExtrapol(){}

    //Accessors
    int                                             diskId() const;
    int                                        trackNumber() const {return _trackNumber;}
    double                                            time() const;
    double                                         timeErr() const;
    double                              pathLengthEntrance() const;
    double                           pathLenghtEntranceErr() const;
    double                                  pathLengthExit() const;
    double                                  fitConsistency() const;
    double                                              t0() const;
    double                                           t0Err() const;
    double                                         tOrigin() const;
    double                                      tOriginErr() const;
    CLHEP::Hep3Vector                     entranceMomentum() const;
    BbrVectorErr                       entranceMomentumErr() const;
    //HepPoint                              entrancePosition() const; #TODO: Check if can delete
    BbrPointErr                        entrancePositionErr() const;
    //HepPoint                                  exitPosition() const;
    BbrPointErr                            exitPositionErr() const;
    CLHEP::Hep3Vector                             momentum() const;
    BbrVectorErr                               momentumErr() const;
    KalRepPtr const&                                   trk() const{ return _trk; }

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

 typedef std::vector<mu2e::TrkToCaloExtrapol> TrkToCaloExtrapolCollection;

} // namespace mu2e

#endif /* TrackCaloMatching_TrkToCaloExtrapol_hh */
