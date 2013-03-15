//
// Create a G4FieldManager object.
//
// $Id: FieldMgr.cc,v 1.4 2013/03/15 19:03:20 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 19:03:20 $
//
// Original author Rob Kutschke
//

// G4 includes.
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

// Mu2e includes
#include "Mu2eG4/inc/FieldMgr.hh"
#include "Mu2eG4/inc/DSGradientField.hh"

using namespace std;

namespace mu2e {

  // Factory method to construct a manager for a uniform magnetic field. See notes in header file.
  std::unique_ptr<FieldMgr> FieldMgr::forUniformField(const G4ThreeVector& fieldValue,
                                                      const G4ThreeVector& mu2eOrigin,
                                                      double stepMinimum ){

    unique_ptr<FieldMgr> mgr(new FieldMgr() );

    mgr->_field       = std::unique_ptr<G4MagneticField>        (new G4UniformMagField   ( fieldValue ));
    mgr->_rhs         = std::unique_ptr<G4Mag_UsualEqRhs>       (new G4Mag_UsualEqRhs    ( mgr->field()) );
    mgr->_integrator  = std::unique_ptr<G4MagIntegratorStepper> (new G4ExactHelixStepper ( mgr->rhs()) );
    mgr->_chordFinder = std::unique_ptr<G4ChordFinder>          (new G4ChordFinder       ( mgr->field(),
                                                                                           stepMinimum,
                                                                                           mgr->integrator()) );
    mgr->_manager     = std::unique_ptr<G4FieldManager>         (new G4FieldManager      ( mgr->field(),
                                                                                           mgr->chordFinder(),
                                                                                           true ));
    return std::move(mgr);
  }

  // Factory method to construct a manager for a gradient magnetic field (in DS3).
  std::unique_ptr<FieldMgr> FieldMgr::forGradientField(double fieldValue,
                                                       double gradient,
                                                       const G4ThreeVector& fieldOrigin,
                                                       double stepMinimum ){

    unique_ptr<FieldMgr> mgr(new FieldMgr() );

    mgr->_field = std::unique_ptr<G4MagneticField>(new DSGradientField( "DS3Grad",
                                                                        fieldOrigin,
                                                                        gradient,
                                                                        fieldValue
                                                                        )
                                                   );
    mgr->_rhs         = std::unique_ptr<G4Mag_UsualEqRhs>       (new G4Mag_UsualEqRhs    ( mgr->field()) );
    mgr->_integrator  = std::unique_ptr<G4MagIntegratorStepper> (new G4ExactHelixStepper ( mgr->rhs()) );
    mgr->_chordFinder = std::unique_ptr<G4ChordFinder>          (new G4ChordFinder       ( mgr->field(),
                                                                                           stepMinimum,
                                                                                           mgr->integrator()) );
    mgr->_manager     = std::unique_ptr<G4FieldManager>         (new G4FieldManager      ( mgr->field(),
                                                                                           mgr->chordFinder(),
                                                                                           true ));
    return std::move(mgr);
  }

  // Release all of the objects that this class owns.
  void FieldMgr::release(){
    _field.release();
    _rhs.release();
    _integrator.release();
    _chordFinder.release();
    _manager.release();
  }

  // Private.  Used only by the factory methods.
  FieldMgr::FieldMgr():
    _field(nullptr),
    _rhs(nullptr),
    _integrator(nullptr),
    _chordFinder(nullptr),
    _manager(nullptr){
  }

} // end namespace mu2e
