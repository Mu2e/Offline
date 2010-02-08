#ifndef DSFIELD_HH
#define DSFIELD_HH

#include <memory>
#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "BFwcont.hh"
#include "G4ThreeVector.hh"

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

#endif
