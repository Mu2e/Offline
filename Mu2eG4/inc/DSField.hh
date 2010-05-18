#ifndef DSFIELD_HH
#define DSFIELD_HH
//
// G4 interface to the Detector Solenoid full magnetic field.
//
// $Id: DSField.hh,v 1.3 2010/05/18 22:33:49 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:33:49 $
//
// Original author Julie Managan and Bob Bernstein

#include <memory>
#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "BFwcont.hh"
#include "G4ThreeVector.hh"

namespace mu2e {

  class DSField: public G4MagneticField {

  public:

    DSField( std::string filename, G4ThreeVector const& origin, 
             int const nx, int const ny, int const nz):
      G4MagneticField(),
      //_p(new BFBMvec(filename,origin,nx,ny,nz)){ 
      _p(new BFwcont(filename,origin,nx,ny,nz)){
    }

    virtual ~DSField(){};
    virtual void GetFieldValue(const G4double Point[4],
                               G4double *Bfield) const;

  private:

    //std::auto_ptr<BFBMvec> _p;
    std::auto_ptr<BFwcont> _p;

  };
}

#endif
