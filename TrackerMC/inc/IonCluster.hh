#ifndef TrackerMC_IonCluster_hh
#define TrackerMC_IonCluster_hh
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"
namespace mu2e {
  struct IonCluster {  // ion charge cluster before drift or amplification
    CLHEP::Hep3Vector _pos; // position of this cluster
    Float_t _charge; // charge of this cluster, in pC.  Note: this is pre-gain!!!
    Float_t _eion; // ionization energy of this cluster, in MeV
    UInt_t _ne; // number of electrons in this cluster
    IonCluster(CLHEP::Hep3Vector const& pos, double charge, double eion, unsigned ne): 
      _pos(pos),_charge(charge),_eion(eion),_ne(ne) {}
    IonCluster() : _charge(0.0), _eion(0.0), _ne(0) {}
  };
}
#endif
