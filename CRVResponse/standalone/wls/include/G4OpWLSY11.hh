/*
GEANT doesn't provide a process which allows multiple absorption and emission spectra, and also uses a quantum yield.
Therefore, two process was developed based of G4OpWLS, one for polystyrene+PPO, and another one for POPOP.
They use different property table names.
The same mechanism was used for Y11 (fibers)
*/


#ifndef G4OpWLSY11_h
#define G4OpWLSY11_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4VWLSTimeGeneratorProfile.hh"

// Class Description:
// Discrete Process -- Bulk absorption of Optical Photons.
// Class inherits publicly from G4VDiscreteProcess
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4VWLSTimeGeneratorProfile;

class G4OpWLSY11 : public G4VDiscreteProcess 
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4OpWLSY11(const G4String& processName = "OpWLSY11",
                         G4ProcessType type = fOptical);
        ~G4OpWLSY11();

private:

        G4OpWLSY11(const G4OpWLSY11 &right);

        //////////////
        // Operators
        //////////////

        G4OpWLSY11& operator=(const G4OpWLSY11 &right);

public:

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

        void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
        // Build the WLS integral table at the right time

        G4double GetMeanFreePath(const G4Track& aTrack,
                                 G4double ,
                                 G4ForceCondition* );
        // Returns the absorption length for bulk absorption of optical
        // photons in media with a specified attenuation length.

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);
        // This is the method implementing bulk absorption of optical
        // photons.

        G4PhysicsTable* GetIntegralTable() const;
        // Returns the address of the WLS integral table.

        void DumpPhysicsTable() const;
        // Prints the WLS integral table.

        void UseTimeProfile(const G4String name);
        // Selects the time profile generator

protected:

        G4VWLSTimeGeneratorProfile* WLSTimeGeneratorProfile;
        G4PhysicsTable* theIntegralTable;

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4OpWLSY11::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
G4PhysicsTable* G4OpWLSY11::GetIntegralTable() const
{
  return theIntegralTable;
}

inline
void G4OpWLSY11::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsOrderedFreeVector *v;
 
  for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
    {
      v = (G4PhysicsOrderedFreeVector*)(*theIntegralTable)[i];
      v->DumpValues();
    }
}

#endif /* G4OpWLSY11_h */
