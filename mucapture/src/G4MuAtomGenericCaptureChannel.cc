// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------
// $Id: G4MuAtomGenericCaptureChannel.cc,v 1.3 2013/10/23 20:50:02 genser Exp $
// $Author: genser $
// $Date: 2013/10/23 20:50:02 $

#include "G4MuAtomGenericCaptureChannel.hh"

#include <cmath>

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

#include "G4MuAtom.hh"
#include "G4Nucleon.hh"
#include "G4Proton.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4NucleiProperties.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4ExcitationHandler.hh"

#include "RandomUtilities.hh"

#include "G4MuAtomTable.hh"
#include "G4MuAtomCaptureRateModel.hh"
#include "G4MuPCaptureChannel.hh"
#include "G4MuDCaptureChannel.hh"
#include "G4MuTCaptureChannel.hh"

G4MuAtomGenericCaptureChannel::
G4MuAtomGenericCaptureChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "generic nuclear capture", 
				  0, verboseLevel) {
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin);
}

G4MuAtomGenericCaptureChannel::
G4MuAtomGenericCaptureChannel(G4MuAtom const* p, G4double rate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "generic nuclear capture", 
				  rate, verboseLevel) {
}

G4MuAtomGenericCaptureChannel* 
G4MuAtomGenericCaptureChannel::Clone(G4MuAtom const* p){
  return new G4MuAtomGenericCaptureChannel(p,verbose);
}

G4DecayProducts* G4MuAtomGenericCaptureChannel::
CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products = new G4DecayProducts;

  G4int const Z = pParticle->GetDefinition()->GetAtomicNumber();
  G4int const A = pParticle->GetDefinition()->GetAtomicMass();

  // FIXME ... what is this doing?
  // This here is just a sanity check
  if( Z<3){
    if( Z==1 ){ // hydrogens
      if( A==1 ) { // protium
      } else if( A==2 ) { // deuterium
      } else if( A==3 ) { // tritium
      }
    } else if( Z==2 ){ // heliums
      if( A==3 ) {
      } else if( A==4 ){
      }
    }
    // FIXME ... this is wrong
    return new G4DecayProducts;
  } 
  // everything Z>2

  // set the parent
  products->SetParentParticle(*pParticle);
  // get the muatom
  G4MuAtom const* muatom = static_cast<G4MuAtom const*>(pParticle->GetDefinition());

  // most of this lifted from class G4MuonMinusNuclearCapture

  // figure initial state parameters
  G4double const mumass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4double const binding = muatom->GetKShellEnergy();
  G4double const muEnergy = mumass + binding;
  G4double const muMom = std::sqrt(binding*(binding+2.*mumass));
  G4ThreeVector const vmu = muMom*GetRandomVec();
  G4LorentzVector const aMuMom(vmu, muEnergy);
  G4LorentzVector const momInitial(0.0,0.0,0.0,muatom->GetPDGMass());
  G4LorentzVector momResidual;
  G4double const residualMass = G4NucleiProperties::GetNuclearMass(A,Z-1);
  G4ReactionProduct* aNu = new G4ReactionProduct(); // FIXME leak
  aNu->SetDefinition( G4NeutrinoMu::NeutrinoMu() );

  // pick random proton inside nucleus 
  G4double eEx;
  static G4Fancy3DNucleus theN;
  do {
    theN.Init(A, Z); 
    G4LorentzVector thePMom;
    G4int theProtonCounter = G4int( Z * G4UniformRand() );
    G4int counter = 0;
    theN.StartLoop();

    G4Nucleon* aNucleon;
    while( (aNucleon = theN.GetNextNucleon()) ) {

      if( aNucleon->GetDefinition() == G4Proton::Proton() ) {
	counter++;
	if(counter == theProtonCounter) {
	  thePMom  = aNucleon->GetMomentum();
	  break;
	}
      }
    }
    // Get the nu momentum in the CMS
    G4LorentzVector theCMS = thePMom + aMuMom;
    G4ThreeVector bst = theCMS.boostVector();

    G4double Ecms = theCMS.mag();
    G4double Enu  = 0.5*(Ecms - CLHEP::neutron_mass_c2*CLHEP::neutron_mass_c2/Ecms);
    eEx = 0.0;

    if(Enu > 0.0) {
      // make the nu, and transform to lab;
      G4ThreeVector nu3Mom = Enu*GetRandomVec();
      G4LorentzVector nuMom(nu3Mom, Enu);

      // nu in lab.
      nuMom.boost(bst);
      aNu->SetTotalEnergy( nuMom.e() );
      aNu->SetMomentum( nuMom.vect() );
    
      // make residual
      momResidual = momInitial - nuMom;

      // Call pre-compound on the rest.
      eEx = momResidual.mag();
      if(GetVerboseLevel() > 1)
	G4cout << "G4MuAtomGenericCaptureChannel::DoMuCapture: " 
	       << " Eex(MeV)= " << (eEx-residualMass)/CLHEP::MeV
	       << " Enu(MeV)= "<<aNu->GetTotalEnergy()/CLHEP::MeV
	       <<G4endl;
    }
  } while(eEx <= residualMass);

  //
  // Start Deexcitation
  //
  G4ThreeVector fromBreit = momResidual.boostVector();
  G4LorentzVector fscm(0.0,0.0,0.0, eEx);
  G4Fragment anInitialState;
  anInitialState.SetA(A);
  anInitialState.SetZ(G4double(Z - 1));
  anInitialState.SetNumberOfParticles(2);
  anInitialState.SetNumberOfCharged(0);
  anInitialState.SetNumberOfHoles(1);
  anInitialState.SetMomentum(fscm);
  static G4ExcitationHandler theHandler;
  G4ReactionProductVector *aPreResult = theHandler.BreakItUp(anInitialState);

  G4ReactionProductVector::iterator ires;
  G4double eBal = muatom->GetPDGMass() - aNu->GetTotalEnergy();
  for(ires=aPreResult->begin(); ires!=aPreResult->end(); ires++) {
    G4LorentzVector itV((*ires)->GetTotalEnergy(), (*ires)->GetMomentum());
    itV.boost(fromBreit);
    (*ires)->SetTotalEnergy(itV.t());
    (*ires)->SetMomentum(itV.vect());
    eBal -= itV.t();
  }

  //
  // fill result
  //
  G4DynamicParticle* p;
  p = new G4DynamicParticle(aNu->GetDefinition(), aNu->GetMomentum());
  products->PushProducts(p);
  for(ires=aPreResult->begin(); ires!=aPreResult->end(); ires++) {
    p = new G4DynamicParticle((*ires)->GetDefinition(), (*ires)->GetMomentum());
    products->PushProducts(p);
  }
    
  delete aPreResult;

  return products;
}
