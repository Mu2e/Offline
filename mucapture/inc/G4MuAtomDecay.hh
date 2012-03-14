// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, cloned from G4Decay,
//      20 February 2010 Kevin Lynch


#ifndef G4MuAtomDecay_h
#define G4MuAtomDecay_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4DecayProcessType.hh"

class G4VExtDecayer;

class G4MuAtomDecay : public G4VRestDiscreteProcess 
{
 // Class Description
  //  This class is a decay process

  public:
    //  Constructors 
    G4MuAtomDecay(const G4String& processName ="MuAtomDecay");

    //  Destructor
    virtual ~G4MuAtomDecay();

  private:
    //  copy constructor
      G4MuAtomDecay(const G4MuAtomDecay &right);

    //  Assignment Operation (generated)
      G4MuAtomDecay & operator=(const G4MuAtomDecay &right);

  public: //With Description
     // G4MuAtomDecay Process has both 
     // PostStepDoIt (for decay in flight) 
     //   and 
     // AtRestDoIt (for decay at rest)
  
     virtual G4VParticleChange *PostStepDoIt(
			     const G4Track& aTrack,
                             const G4Step& aStep
                            );

     virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );

     virtual void BuildPhysicsTable(const G4ParticleDefinition&); 
     // In G4MuAtomDecay, thePhysicsTable stores values of
    //    beta * std::sqrt( 1 - beta*beta) 
    //  as a function of normalized kinetic enregy (=Ekin/mass),
    //  becasuse this table is universal for all particle types,


    virtual G4bool IsApplicable(const G4ParticleDefinition&);
    // returns "true" if the decay process can be applied to
    // the particle type. 
 
  protected: // With Description
    virtual G4VParticleChange* DecayIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    );
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

    // Set daughter polarization
    //  NO OPERATION in the base class of G4MuAtomDecay 
    virtual void DaughterPolarization(const G4Track& aTrack,
			      G4DecayProducts* products);

 public:
    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            );

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

  protected: // With Description
    // GetMeanFreePath returns ctau*beta*gamma for decay in flight 
    // GetMeanLifeTime returns ctau for decay at rest
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double   previousStepSize,
                              G4ForceCondition* condition
                             );

    virtual G4double GetMeanLifeTime(const G4Track& aTrack,
                              G4ForceCondition* condition
                            );

   public: //With Description
     virtual void StartTracking(G4Track*);
     virtual void EndTracking();
      // inform Start/End of tracking for each track to the physics process 

   public: //With Description
     void SetExtDecayer(G4VExtDecayer*);
     const G4VExtDecayer* GetExtDecayer() const;
     // Set/Get External Decayer
   
    G4double GetRemainderLifeTime() const;  
    //Get Remainder of life time at rest decay 

  public:
     void  SetVerboseLevel(G4int value);
     G4int GetVerboseLevel() const;

  protected:
     G4int verboseLevel;
     // controle flag for output message
     //  0: Silent
     //  1: Warning message
     //  2: More

  protected:
    // HighestValue.
    const G4double HighestValue;
 
    // Remainder of life time at rest
    G4double                 fRemainderLifeTime;
  
    // ParticleChange for decay process
    G4ParticleChangeForDecay fParticleChangeForDecay;
    
    // External Decayer
    G4VExtDecayer*    pExtDecayer;
};

inline
 void  G4MuAtomDecay::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4MuAtomDecay::GetVerboseLevel() const { return verboseLevel; }

inline  
  G4VParticleChange* G4MuAtomDecay::AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}

inline  
  G4VParticleChange* G4MuAtomDecay::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step&  aStep
			    )
{
  return DecayIt(aTrack, aStep);
}

inline
 const G4VExtDecayer* G4MuAtomDecay::GetExtDecayer() const
{
  return pExtDecayer;
}

inline
 G4double G4MuAtomDecay::GetRemainderLifeTime() const 
{
  return fRemainderLifeTime;
}

#endif










