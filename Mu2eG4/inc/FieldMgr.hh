#ifndef Mu2eG4_FieldMgr_hh
#define Mu2eG4_FieldMgr_hh
//
// Create a G4FieldManager object. Provide accessors to the field manager
// and to the parts from which it is made.
//
// $Id: FieldMgr.hh,v 1.10 2013/10/25 21:42:02 genser Exp $
// $Author: genser $
// $Date: 2013/10/25 21:42:02 $
//
// Original author Rob Kutschke
//
// Notes
// 1) In the general case, construction of a G4FieldManager is 5 step process.
//    This class does this 5 step process in 2 special cases.
//    a) A uniform magnetic field.
//    b) A general magnetic field described by a field map held by the
//       Mu2e GeometeryService.
//
// 2) Each of the 5 steps involves new'ing an object.  Each of these objects
//    must have a lifetime matched to the lifetime of the G4Geometry.
//
// 3) This class also manages the lifetime of the 5 objects: they will have
//    a lifetime equal to the lifetime of an object in this class.
//
// 4) For case of 1b), this class provides the option of supplying different
//    G4StepperIntegrator objects.  Examples include G4ExplicitEuler,
//    G4ImplicitEuler, G4ClassicalRK4, G4CashKarpRKF45 and others.
//
// 5) To implement the options in 5, this class uses templates, rather
//    than a base class and many concrete classes.
//
// 6) For case 1b), because of limitations in C++ syntax, this class does
//    not provide a constructor for FieldMgr objects.  Instead it provides
//    provides a templated factory method to create instances of FieldMgr objects.
//    It is possible to write a constructor but it requires a confusing
//    extra argument in order to make the template type deduction work correctly.
//
// 7) For case 1a), it would be easy to write a constructor but I wrote
//    a factory method to maintain the symmetry of the class.
//
// 8) The magic number for the default of stepMinimum comes from the source
//    code for G4Chordfinder, v1.21
//
#include <memory>
#include <string>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// G4 includes
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"
#include "Mu2eG4/inc/Mu2eGlobalField.hh"

namespace mu2e {

  class FieldMgr{

  public:

    // See private section for c'tors and assignment operator.
    // Accept compiler generated d'tor.

    // Accessors: return non-owning bare pointers to the bits and pieces that make up the field manager.
    // G4 requires bare pointers so that is what we provide.
    G4MagneticField*        field()       { return _field.get(); }
    G4Mag_EqRhs*            rhs()         { return _rhs.get(); }
    G4MagIntegratorStepper* integrator()  { return _integrator.get(); }
    G4ChordFinder*          chordFinder() { return _chordFinder.get(); }
    G4FieldManager*         manager()     { return _manager.get(); }

    // Factory method to construct a manager for a uniform magnetic field.  See Note 8.
    static std::unique_ptr<FieldMgr> forUniformField(const G4ThreeVector& fieldValue,
                                                   const G4ThreeVector& mu2eOrigin,
                                                   double stepMinimum=1.0e-2*CLHEP::mm);

    // Factory method to construct a manager for a gradient magnetic field.
    static std::unique_ptr<FieldMgr> forGradientField(double fieldValue,
						   double gradient,
						   const G4ThreeVector& fieldOrigin,
                                                   double stepMinimum=1.0e-2*CLHEP::mm);

    // Factory method to construct a manager for a magnetic field described by a Mu2e field map
    // and will a user supplied G4IntegratorStepper.  See Note 8.
    template <class INTEGRATOR>
    static std::unique_ptr<FieldMgr> forMappedField(const std::string& fieldName,
                                                  const G4ThreeVector& mu2eOrigin,
                                                  double stepMinimum=1.0e-2*CLHEP::mm){

      std::unique_ptr<FieldMgr> mgr(new FieldMgr() );

      mgr->_field       = std::unique_ptr<G4MagneticField>        (new Mu2eGlobalField  ( mu2eOrigin) );
      mgr->_rhs         = std::unique_ptr<G4Mag_UsualEqRhs>       (new G4Mag_UsualEqRhs ( mgr->field()) );
      mgr->_integrator  = std::unique_ptr<G4MagIntegratorStepper> (new INTEGRATOR       ( mgr->rhs()) );
      mgr->_chordFinder = std::unique_ptr<G4ChordFinder>          (new G4ChordFinder    ( mgr->field(),
                                                                                       stepMinimum,
                                                                                       mgr->integrator()) );
      mgr->_manager     = std::unique_ptr<G4FieldManager>         (new G4FieldManager   ( mgr->field(),
                                                                                       mgr->chordFinder(),
                                                                                       true ));
      return mgr;
    }

    // Release all of the objects that this class owns.
    void release();

  private:

    // The parts that make up a field manager.
    std::unique_ptr<G4MagneticField>        _field;
    std::unique_ptr<G4Mag_EqRhs>            _rhs;
    std::unique_ptr<G4MagIntegratorStepper> _integrator;
    std::unique_ptr<G4ChordFinder>          _chordFinder;
    std::unique_ptr<G4FieldManager>         _manager;

    // Only the factory methods should make new instances of this class.
    FieldMgr();

    // This class should never be copied because the unique_ptr obeys move semantics.
    FieldMgr( const FieldMgr& );
    FieldMgr& operator=(const FieldMgr& );

  };

} //end namespace mu2e

#endif /* Mu2eG4_FieldMgr_hh */
