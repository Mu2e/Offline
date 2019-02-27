/*
Cerenkov process which doesn't assume that the the index of refraction increases with increasing energy.
In this more general case, the integrals, and cut-offs need to be treated differently.
*/

#ifndef G4Cerenkov_h
#define G4Cerenkov_h 1

/////////////
// Includes
/////////////

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

#include <map>

// Class Description:
// Discrete Process -- Generation of Cerenkov Photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4CerenkovNew : public G4VProcess
{

public:

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

	G4CerenkovNew(const G4String& processName = "CerenkovNew", 
                            G4ProcessType type = fElectromagnetic);
	~G4CerenkovNew();

        G4CerenkovNew(const G4CerenkovNew &right);

private:

        //////////////
        // Operators
        //////////////

        G4CerenkovNew& operator=(const G4CerenkovNew &right);

public:

        ////////////
        // Methods
        ////////////

        static G4CerenkovNew* Instance() {return _fgInstance;}

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable', for all charged particles
        // except short-lived particles.

        void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
        // Build table at a right time

        G4double GetMeanFreePath(const G4Track& aTrack,
                                 G4double ,
                                 G4ForceCondition* );
        // Returns the discrete step limit and sets the 'StronglyForced'
        // condition for the DoIt to be invoked at every step.

        G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                      G4double ,
                                                      G4ForceCondition* );
        // Returns the discrete step limit and sets the 'StronglyForced'
        // condition for the DoIt to be invoked at every step.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
					const G4Step&  aStep);
        // This is the method implementing the Cerenkov process.

        //  no operation in  AtRestDoIt and  AlongStepDoIt
        virtual G4double AlongStepGetPhysicalInteractionLength(
                               const G4Track&,
                               G4double  ,
                               G4double  ,
                               G4double& ,
                               G4GPILSelection*
                              ) { return -1.0; };

        virtual G4double AtRestGetPhysicalInteractionLength(
                               const G4Track& ,
                               G4ForceCondition*
                              ) { return -1.0; };

        //  no operation in  AtRestDoIt and  AlongStepDoIt
        virtual G4VParticleChange* AtRestDoIt(
                               const G4Track& ,
                               const G4Step&
                              ) {return 0;};

        virtual G4VParticleChange* AlongStepDoIt(
                               const G4Track& ,
                               const G4Step&
                              ) {return 0;};

        void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any 
        // produced Cerenkov photons are tracked next. When all have 
        // been tracked, the tracking of the primary resumes.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.

        void SetMaxBetaChangePerStep(const G4double d);
        // Set the maximum allowed change in beta = v/c in % (perCent)
        // per step.

        G4double GetMaxBetaChangePerStep() const;
        // Returns the maximum allowed change in beta = v/c in % (perCent)

        void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
        // Set the maximum number of Cerenkov photons allowed to be 
        // generated during a tracking step. This is an average ONLY; 
        // the actual number will vary around this average. If invoked, 
        // the maximum photon stack will roughly be of the size set.
        // If not called, the step is not limited by the number of 
        // photons generated.

        G4int GetMaxNumPhotonsPerStep() const;
        // Returns the maximum number of Cerenkov photons allowed to be
        // generated during a tracking step.

        G4PhysicsTable* GetPhysicsTable() const;
        // Returns the address of the physics table.

        void DumpPhysicsTable() const;
        // Prints the physics table.

	G4double GetAverageNumberOfPhotons(const G4double charge,
                                const G4double beta,
		    		const G4Material *aMaterial,
				G4MaterialPropertyVector* Rindex) const;

private:

        void BuildThePhysicsTable();

	/////////////////////
	// Helper Functions
	/////////////////////

        void GetMinMaxRindex(G4MaterialPropertyVector* Rindex, double &RImin, double &RImax) const;

        void GetEnergyIntervals(G4MaterialPropertyVector* Rindex, double BetaInverse, 
                                std::vector<std::pair<double,double> > &energyIntervals) const;

        ///////////////////////
        // Class Data Members
        ///////////////////////

protected:

        G4PhysicsTable* thePhysicsTable;
        //  A Physics Table can be either a cross-sections table or
        //  an energy table (or can be used for other specific
        //  purposes).

        std::map<int,double> minRindex;
        std::map<int,double> maxRindex;

private:

        static G4CerenkovNew*  _fgInstance;  

        G4bool fTrackSecondariesFirst;
        G4double fMaxBetaChange;
        G4int  fMaxPhotons;
};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4CerenkovNew::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline
G4double G4CerenkovNew::GetMaxBetaChangePerStep() const
{
        return fMaxBetaChange;
}

inline
G4int G4CerenkovNew::GetMaxNumPhotonsPerStep() const
{
        return fMaxPhotons;
}

inline
void G4CerenkovNew::DumpPhysicsTable() const
{
        G4int PhysicsTableSize = thePhysicsTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
        {
        	v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
        	v->DumpValues();
        }
}

inline
G4PhysicsTable* G4CerenkovNew::GetPhysicsTable() const
{
  return thePhysicsTable;
}

#endif /* G4CerenkovNew_h */
