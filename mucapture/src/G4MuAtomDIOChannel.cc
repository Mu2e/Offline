// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuAtomDIOChannel.hh"
#include "G4MuAtom.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "RandomUtilities.hh"

// FIXME ... I don't like this implementation.  This stuff shouldn't be
// in the constructor. 
G4MuAtomDIOChannel::
G4MuAtomDIOChannel(G4String const& parentName, G4double BR) : 
  G4VMuAtomDecayChannel("MuAtom DIO", 1)
{
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *def = table->FindParticle(parentName);
  if( def->GetParticleType() == "MuAtom" ) {
    SetBR(BR);
    SetParent(parentName);
  } else {
    if( GetVerboseLevel()>0 ){
      G4cout << "G4MuAtomDIOChannel::constructor : ";
      G4cout << " parent particle is not a muAtom, but ";
      G4cout << parentName << '\n';
    }
  }
  SetNumberOfDaughters(3);
  SetDaughter(0,"e-");
  SetDaughter(1,"anti_nu_e");
  SetDaughter(2,"nu_mu");
}

G4MuAtomDIOChannel::~G4MuAtomDIOChannel(){
}

G4MuAtomDIOChannel* G4MuAtomDIOChannel::Clone(){
  G4MuAtomDIOChannel *chan = new G4MuAtomDIOChannel(*this);
  return chan;
}

// FIXME ... check the physics carefully ... I'm not sure it's right.

// partially snarfed from G4MuonDecayChannel, and augmented by
// G4MuMinusCaptureCascade::DoBoundMuonDecay 
G4DecayProducts* G4MuAtomDIOChannel::DecayIt(G4double){

  // G4MuonDecayChannel::DecayIt
  if( parent == 0 )
    FillParent();
  if( daughters == 0 )
    FillDaughters();
  
  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }
  
  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  G4MuAtom* muatom = static_cast<G4MuAtom*>(parent);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double const Emass = G4Electron::Definition()->GetPDGMass();
  G4double const MuMass = G4MuonMinus::Definition()->GetPDGMass();

  // G4MuMinusCaptureCascade::DoBoundMuonDecay
  // Simulation on Decay of mu- on a K-shell of the muonic atom
  G4double xmax = ( 1.0 + Emass*Emass/ (MuMass*MuMass) );
  G4double xmin = 2.0*Emass/MuMass;
  G4double KEnergy = muatom->GetKShellEnergy();

  G4double pmu = std::sqrt(KEnergy*(KEnergy + 2.0*MuMass));
  G4double emu = KEnergy + MuMass;
  G4ThreeVector moment = GetRandomVec();
  G4LorentzVector MU(pmu*moment,emu);
  G4ThreeVector bst = MU.boostVector();

  // Now, do the decay in the muon rest frame ... we'll boost to the
  // moving frame after

  G4double Eelect, Pelect, x, ecm;
  G4LorentzVector EL, NN, N1, N2;
  // Calculate electron energy ... FIXME ... I changed this from
  // G4MuMinusCaptureCascade, and I'm not sure I'm getting the right
  // results.  Check it out...
  do {
    do {
      x = xmin + (xmax-xmin)*G4UniformRand();
    } while (G4UniformRand() > (3.0 - 2.0*x)*x*x );
    Eelect = x*MuMass*0.5;
    Pelect = 0.0;
    if(Eelect > Emass) { 
      Pelect = std::sqrt( Eelect*Eelect - Emass*Emass );
    } else {
      Pelect = 0.0;
      Eelect = Emass;
    }
    G4ThreeVector e_mom = GetRandomVec();
    EL = G4LorentzVector(Pelect*e_mom,Eelect);
    //    EL.boost(bst);
    //    Eelect = EL.e() - Emass - 2.0*KEnergy;
    //
    // Calculate rest frame parameters of 2 neutrinos
    //
    //    NN = MU - EL;
    //    ecm = NN.mag2();
  } while (Eelect < 0.0);// || ecm < 0.0);

  NN = G4LorentzVector(G4ThreeVector(),MuMass) - EL;
  ecm = NN.mag2();
  ecm = 0.5*std::sqrt(ecm);
  G4ThreeVector p1 = ecm * GetRandomVec();
  N1 = G4LorentzVector(p1,ecm);
  N2 = NN - N1;

  // create dynamic electron
  EL.boost(bst);
  G4DynamicParticle* dynElec = new G4DynamicParticle(daughters[0], EL.vect());
  products->PushProducts(dynElec);

  // create dynamic neutrinoes
  N1.boost(bst);
  G4DynamicParticle* dynAnti = new G4DynamicParticle(daughters[1], N1.vect());
  products->PushProducts(dynAnti);

  N2.boost(bst);
  G4DynamicParticle* dynNeut = new G4DynamicParticle(daughters[2], N2.vect());
  products->PushProducts(dynNeut);

  // FIXME! where's the nucleus?  that should be here too, especially
  // for DIO for a system in motion
  
  return products;
}


