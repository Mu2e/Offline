//
// Class Description:
//
// G4 class that takes produces a conversion electron from an at rest
// muon. Configurable parameters include endpoint of momentum, limits of polar
// and azimuthal angle.
//
// $Id: MuonMinusConversionAtRest.hh,v 1.4 2011/05/18 05:10:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 05:10:17 $
//
//-----------------------------------------------------------------------------

#ifndef muMinusConversionAtRest_h
#define muMinusConversionRest_h 1

#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"

// Mu2e includes
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandomEngine.h"

namespace fhicl { class ParameterSet; }

namespace mu2e {
   class SimpleConfig;

class muMinusConversionAtRest : public G4VRestProcess

{
public:

  muMinusConversionAtRest( const SimpleConfig& config, const G4String& processName ="muMinusConversionAtRest", G4ProcessType   aType = fHadronic );

  muMinusConversionAtRest( const fhicl::ParameterSet& config, const G4String& processName ="muMinusConversionAtRest", G4ProcessType   aType = fHadronic );

   ~muMinusConversionAtRest();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);

  G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*)
  {return 0;};

private:

  // hide assignment operator as private
  muMinusConversionAtRest& operator=(const muMinusConversionAtRest &right);
  muMinusConversionAtRest(const muMinusConversionAtRest& );

  // Conversion momentum.
  G4double _p;

  // Limits on the generated direction.
  G4double _czmin;
  G4double _czmax;
  G4double _phimin;
  G4double _phimax;

  //Utility to generate direction of the momentum
  mu2e::RandomUnitSphere    _randomUnitSphere;

  G4bool isInitialised;

};


}

#endif






