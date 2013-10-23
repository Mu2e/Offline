// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, cloned from G4Decay
//      20 February 2010 Kevin Lynch

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

#include "G4MuAtom.hh"
#include "G4MuAtomDecay.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4MuAtomDecayTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4VExtDecayer.hh"

// constructor
G4MuAtomDecay::G4MuAtomDecay(const G4String& processName)
                               :G4VRestDiscreteProcess(processName, fDecay),
				verboseLevel(1),
                                HighestValue(20.0),
                                pExtDecayer(0)
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(DECAY));

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4MuAtomDecay  constructor " << "  Name:" << processName << G4endl;
  }
#endif

  pParticleChange = &fParticleChangeForDecay;
}

G4MuAtomDecay::~G4MuAtomDecay()
{
  if (pExtDecayer) {
    delete pExtDecayer;
  }
}

G4bool G4MuAtomDecay::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   // check if the particle is stable?
   if (aParticleType.GetPDGLifeTime() <0.0) {
     return false;
   } else if (aParticleType.GetPDGMass() <= 0.0*CLHEP::MeV) {
     return false;
   } else if (aParticleType.GetParticleType() == "MuAtom") {
     return true; 
   } else {
     return false;
   }
}

G4double G4MuAtomDecay::GetMeanLifeTime(const G4Track& aTrack  ,
                                  G4ForceCondition*)
{
   // returns the mean free path in GEANT4 internal units
   G4double meanlife;

   // get particle 
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aLife = aParticleDef->GetPDGLifeTime();

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
     meanlife = DBL_MAX;
    
   } else {
     meanlife = aLife;
   }

#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
     G4cout << "mean life time: "<< meanlife/CLHEP::ns << "[ns]" << G4endl;
   }
#endif

   return  meanlife;
}

G4double G4MuAtomDecay::GetMeanFreePath(const G4Track& aTrack,G4double, G4ForceCondition*)
{
   // get particle 
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();
   G4double aMass = aParticle->GetMass();
   G4double aLife = aParticleDef->GetPDGLifeTime();


    // returns the mean free path in GEANT4 internal units
   G4double pathlength;
   G4double aCtau = CLHEP::c_light * aLife;

   // check if the particle is stable?
   if (aParticleDef->GetPDGStable()) {
     pathlength = DBL_MAX;

   //check if the particle has very short life time ?
   } else if (aCtau < DBL_MIN) { 
     pathlength =  DBL_MIN;
 
   } else {
    //calculate the mean free path 
    // by using normalized kinetic energy (= Ekin/mass)
     G4double   rKineticEnergy = aParticle->GetKineticEnergy()/aMass; 
     if ( rKineticEnergy > HighestValue) {
       // gamma >>  1
       pathlength = ( rKineticEnergy + 1.0)* aCtau;
     } else if ( rKineticEnergy < DBL_MIN ) {
       // too slow particle
#ifdef G4VERBOSE
       if (GetVerboseLevel()>1) {
	 G4cout << "G4MuAtomDecay::GetMeanFreePath()   !!particle stops!!";
         G4cout << aParticleDef->GetParticleName() << G4endl;
	 G4cout << "KineticEnergy:" << aParticle->GetKineticEnergy()/CLHEP::GeV <<"[GeV]";
       }
#endif
       pathlength = DBL_MIN;
     } else {
       // beta <1 
       pathlength = (aParticle->GetTotalMomentum())/aMass*aCtau ;
     }
   }
  return  pathlength;
}

void G4MuAtomDecay::BuildPhysicsTable(const G4ParticleDefinition&)
{
  return;
}

G4VParticleChange* G4MuAtomDecay::DecayIt(const G4Track& aTrack, const G4Step& )
{
  // The DecayIt() method returns by pointer a particle-change object.
  // Units are expressed in GEANT4 internal units.

  //   Initialize ParticleChange
  //     all members of G4VParticleChange are set to equal to 
  //     corresponding member in G4Track
  fParticleChangeForDecay.Initialize(aTrack);

  // get particle 
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4MuAtom* aParticleDef = 
    static_cast<G4MuAtom*>(aParticle->GetDefinition());

  // check if  the particle is stable
  if (aParticleDef->GetPDGStable()) return &fParticleChangeForDecay ;
 

  //check if thePreAssignedDecayProducts exists
  const G4DecayProducts* o_products = (aParticle->GetPreAssignedDecayProducts());
  G4bool isPreAssigned = (o_products != 0);   
  G4DecayProducts* products = 0;

  // decay table
  G4MuAtomDecayTable   *decaytable = aParticleDef->GetMuAtomDecayTable();
 
  // check if external decayer exists
  G4bool isExtDecayer = (decaytable == 0) && (pExtDecayer !=0);

  // Error due to NO Decay Table 
  if ( (decaytable == 0) && !isExtDecayer &&!isPreAssigned ){
    if (GetVerboseLevel()>0) {
      G4cout <<  "G4MuAtomDecay::DoIt  : decay table not defined  for ";
      G4cout << aParticle->GetDefinition()->GetParticleName()<< G4endl;
    }
    G4Exception( "G4MuAtomDecay::DecayIt ",
                 "MuAtomDecay0001",JustWarning, 
                 "Decay table is not defined");

    fParticleChangeForDecay.SetNumberOfSecondaries(0);
    // Kill the parent particle
    fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
    fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.0); 
    
    ClearNumberOfInteractionLengthLeft();
    return &fParticleChangeForDecay ;
  }

  if (isPreAssigned) {
    // copy decay products 
    products = new G4DecayProducts(*o_products); 
  } else if ( isExtDecayer ) {
    // decay according to external decayer
    products = pExtDecayer->ImportDecayProducts(aTrack);
  } else {
    // decay acoording to decay table
    // choose a decay channel
    G4VMuAtomDecayChannel *decaychannel = decaytable->SelectADecayChannel();
    if (decaychannel == 0 ){
      // decay channel not found
      //      G4Exception("G4MuAtomDecay::DoIt  : can not determine decay channel ");
      G4Exception( "G4MuAtomDecay::DoIt",
                   "MuAtomDecay0002",FatalException,
                   "can not determine decay channel");
    } else {
      // execute DecayIt() 
#ifdef G4VERBOSE
      G4int temp = decaychannel->GetVerboseLevel();
      if (GetVerboseLevel()>1) {
	G4cout << "G4MuAtomDecay::DoIt  : selected decay channel  addr:" << decaychannel <<G4endl;
	decaychannel->SetVerboseLevel(GetVerboseLevel());
      }
#endif
      products = decaychannel->DecayIt(aParticle->GetMass());
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
	decaychannel->SetVerboseLevel(temp);
      }
#endif
#ifdef G4VERBOSE
      if (GetVerboseLevel()>2) {
	if (! products->IsChecked() ) products->DumpInfo();
      }
#endif
    }
  }
  
  // get parent particle information ...................................
  G4double   ParentEnergy  = aParticle->GetTotalEnergy();
  G4double   ParentMass    = aParticle->GetMass();
  if (ParentEnergy < ParentMass) {
    ParentEnergy = ParentMass;
    if (GetVerboseLevel()>0) {
      G4cout << "G4MuAtomDecay::DoIt  : Total Energy is less than its mass" << G4endl;
      G4cout << " Particle: " << aParticle->GetDefinition()->GetParticleName();
      G4cout << " Energy:"    << ParentEnergy/CLHEP::MeV << "[MeV]";
      G4cout << " Mass:"    << ParentMass/CLHEP::MeV << "[MeV]";
      G4cout << G4endl;
    }
  }

  G4ThreeVector ParentDirection(aParticle->GetMomentumDirection());

  //boost all decay products to laboratory frame
  G4double energyDeposit = 0.0;
  G4double finalGlobalTime = aTrack.GetGlobalTime();
  if (aTrack.GetTrackStatus() == fStopButAlive ){
    // AtRest case
    finalGlobalTime += fRemainderLifeTime;
    energyDeposit += aParticle->GetKineticEnergy();
    if (isPreAssigned) products->Boost( ParentEnergy, ParentDirection);
  } else {
    // PostStep case
    if (!isExtDecayer) products->Boost( ParentEnergy, ParentDirection);
  }

   // set polarization for daughter particles
   DaughterPolarization(aTrack, products);


  //add products in fParticleChangeForDecay
  G4int numberOfSecondaries = products->entries();
  fParticleChangeForDecay.SetNumberOfSecondaries(numberOfSecondaries);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4MuAtomDecay::DoIt  : Decay vertex :";
    G4cout << " Time: " << finalGlobalTime/CLHEP::ns << "[ns]";
    G4cout << " X:" << (aTrack.GetPosition()).x() /CLHEP::cm << "[cm]";
    G4cout << " Y:" << (aTrack.GetPosition()).y() /CLHEP::cm << "[cm]";
    G4cout << " Z:" << (aTrack.GetPosition()).z() /CLHEP::cm << "[cm]";
    G4cout << G4endl;
    G4cout << "G4MuAtomDecay::DoIt  : decay products in Lab. Frame" << G4endl;
    products->DumpInfo();
  }
#endif
  G4int index;
  G4ThreeVector currentPosition;
  const G4TouchableHandle thand = aTrack.GetTouchableHandle();
  for (index=0; index < numberOfSecondaries; index++)
  {
     // get current position of the track
     currentPosition = aTrack.GetPosition();
     // create a new track object
     G4Track* secondary = new G4Track( products->PopProducts(),
				      finalGlobalTime ,
				      currentPosition );
     // switch on good for tracking flag
     secondary->SetGoodForTrackingFlag();
     secondary->SetTouchableHandle(thand);
     // add the secondary track in the List
     fParticleChangeForDecay.AddSecondary(secondary);
  }
  delete products;

  // Kill the parent particle
  fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(energyDeposit); 
  fParticleChangeForDecay.ProposeGlobalTime( finalGlobalTime );
  // Clear NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  return &fParticleChangeForDecay ;
} 

void G4MuAtomDecay::DaughterPolarization(const G4Track& , G4DecayProducts* )
{
}



void G4MuAtomDecay::StartTracking(G4Track*)
{
  currentInteractionLength = -1.0;
  ResetNumberOfInteractionLengthLeft();
 
  fRemainderLifeTime = -1.0;
}

void G4MuAtomDecay::EndTracking()
{
  // Clear NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  currentInteractionLength = -1.0;
}


G4double G4MuAtomDecay::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
 
   // condition is set to "Not Forced"
  *condition = NotForced;

   // pre-assigned Decay time
  G4double pTime = track.GetDynamicParticle()->GetPreAssignedDecayProperTime();
  G4double aLife = track.GetDynamicParticle()->GetDefinition()->GetPDGLifeTime();

  if (pTime < 0.) {
    // normal case 
    if ( previousStepSize > 0.0){
      // subtract NumberOfInteractionLengthLeft 
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft<0.){
	theNumberOfInteractionLengthLeft=CLHEP::perMillion;
      }
      fRemainderLifeTime = theNumberOfInteractionLengthLeft*aLife;
    }
    // get mean free path
    currentInteractionLength = GetMeanFreePath(track, previousStepSize, condition);
    
#ifdef G4VERBOSE
    if ((currentInteractionLength <=0.0) || (verboseLevel>2)){
      G4cout << "G4MuAtomDecay::PostStepGetPhysicalInteractionLength " << G4endl;
      track.GetDynamicParticle()->DumpInfo();
      G4cout << " in Material  " << track.GetMaterial()->GetName() <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/CLHEP::cm << "[cm]" <<G4endl;
    }
#endif

    G4double value;
    if (currentInteractionLength <DBL_MAX) {
      value = theNumberOfInteractionLengthLeft * currentInteractionLength;
    } else {
      value = DBL_MAX;
    }

    return value;

  } else {
    //pre-assigned Decay time case
    // reminder proper time
    fRemainderLifeTime = pTime - track.GetProperTime();
    if (fRemainderLifeTime <= 0.0) fRemainderLifeTime = DBL_MIN;
    
    G4double  rvalue=0.0; 
    // use pre-assigned Decay time to determine PIL
    if (aLife>0.0) {
      // ordinary particle
      rvalue = (fRemainderLifeTime/aLife)*GetMeanFreePath(track, previousStepSize, condition);
    } else {
     // shortlived particle
      rvalue = CLHEP::c_light * fRemainderLifeTime;
     // by using normalized kinetic energy (= Ekin/mass)
     G4double   aMass =  track.GetDynamicParticle()->GetMass();
     rvalue *= track.GetDynamicParticle()->GetTotalMomentum()/aMass;
    }
    return rvalue;
  }
}

G4double G4MuAtomDecay::AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            )
{
     // condition is set to "Not Forced"
  *condition = NotForced;

  G4double pTime = track.GetDynamicParticle()->GetPreAssignedDecayProperTime();
  if (pTime >= 0.) {
    fRemainderLifeTime = pTime - track.GetProperTime();
    if (fRemainderLifeTime <= 0.0) fRemainderLifeTime = DBL_MIN;
  } else {
    fRemainderLifeTime = 
      theNumberOfInteractionLengthLeft * GetMeanLifeTime(track, condition);
  }
  return fRemainderLifeTime;
}


void G4MuAtomDecay::SetExtDecayer(G4VExtDecayer* val)
{
  pExtDecayer = val;

  // set Process Sub Type
  if ( pExtDecayer !=0 ) {
    SetProcessSubType(static_cast<int>(DECAY_External));
  }
}
