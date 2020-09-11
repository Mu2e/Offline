//
// ProtonAbsorber description for TrkExt
//
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtProtonAbsorber_HH
#define TrkExtProtonAbsorber_HH

#include <string>
#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtShape.hh"
#include "TrkExt/inc/TrkExtMaterial.hh"

namespace mu2e {

  class TrkExtProtonAbsorber : public TrkExtShape, public TrkExtMaterial 
  {

  public:
    TrkExtProtonAbsorber() ;
    ~TrkExtProtonAbsorber() { }

    void initialize () ;
    bool contains (CLHEP::Hep3Vector& p) ;

  private:
    double z0, z1;
    double r0in, r0out, r1in, r1out;
    double slope_in, slope_out;
    bool valid;
    std::string name;

  };



} // end namespace mu2e


#endif
