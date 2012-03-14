// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 30 2010
// ----------------------------------------------------------------

#include "G4DMuDFusionChannels.hh"

#include "G4MuAtom.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Triton.hh"

#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4DecayProducts.hh"

#include <sstream>
#include <strstream>

// D_mu_D -> He3 + n + mu

G4DMuDFusionHe3Channel::G4DMuDFusionHe3Channel(G4MuMolecule const* p, G4double fusionRate, 
					       G4int verboseLevel) :
  G4VMuMoleculeCaptureKineticsChannel(p, "D_mu_D -> He3 + n + mu fusion", 
				      fusionRate, verboseLevel),
  psdc(p->GetParticleName(), 1., 2, "neutron", "He3"){
  CheckIsApplicable();
}


G4DecayProducts* G4DMuDFusionHe3Channel::CaptureIt(G4DynamicParticle const* p){
  // FIXME ... not handling the molecular binding energy correctly
  G4MuonMinus *muon = G4MuonMinus::Definition();
  // first, the fusion products ...
  G4DecayProducts *products = 
    psdc.DecayIt( p->GetMass() - muon->GetPDGMass() );
  // ... then, the spectator muon ...
  products->PushProducts( new G4DynamicParticle(muon, G4ThreeVector()) );
  // ... finally, the parent
  products->SetParentParticle( *p );
  return products;
}


void G4DMuDFusionHe3Channel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "D_mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o << "The channel " << GetChannelName() 
      << " is only applicable to D_mu_D molecules, not "
      << part->GetParticleName();
    G4Exception("G4DMuDFusionHe3Channel::CheckIsApplicable",
                "DMUDF0001", FatalException, 
                o.str());
  }
}


// D_mu_D -> mu_He3 + n

G4DMuDFusionMuHe3Channel::
G4DMuDFusionMuHe3Channel(G4MuMolecule const* p, G4double fusionRate, 
			 G4int verboseLevel) :
  G4VMuMoleculeCaptureKineticsChannel(p, "D_mu_D -> mu_He3 + n fusion", 
				      fusionRate, verboseLevel),
  psdc(p->GetParticleName(), 1., 2, "neutron"){
  CheckIsApplicable();
  // mu_He3
  psdc.SetDaughter(1,G4MuAtom::Definition(2,3));
}


G4DecayProducts* G4DMuDFusionMuHe3Channel::CaptureIt(G4DynamicParticle const* p){
  // FIXME ... the sticking case needs the muon mass in the parent
  // particle.  The binding energy is still not treated properly...
  G4DecayProducts *products = psdc.DecayIt( p->GetMass() );
  // ... then, the parent
  products->SetParentParticle( *p );
  return products;
}


void G4DMuDFusionMuHe3Channel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "D_mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o << "The channel " << GetChannelName() 
      << " is only applicable to D_mu_D molecules, not "
      << part->GetParticleName();
    G4Exception("G4DMuDFusionMuHe3Channel::CheckIsApplicable",
                "DMUDF0002", FatalException, 
                o.str());
  }
}


// D_mu_D -> T + p + mu

G4DMuDFusionTChannel::G4DMuDFusionTChannel(G4MuMolecule const* p, G4double fusionRate, 
					       G4int verboseLevel) :
  G4VMuMoleculeCaptureKineticsChannel(p, "D_mu_D -> T + n + mu fusion", 
				      fusionRate, verboseLevel),
  psdc(p->GetParticleName(), 1., 2, "proton", "triton"){
  CheckIsApplicable();
}


G4DecayProducts* G4DMuDFusionTChannel::CaptureIt(G4DynamicParticle const* p){
  // FIXME ... not handling the molecular binding energy correctly
  G4MuonMinus *muon = G4MuonMinus::Definition();
  // first, the fusion products ...
  G4DecayProducts *products = 
    psdc.DecayIt( p->GetMass() - muon->GetPDGMass() );
  // ... then, the spectator muon ...
  products->PushProducts( new G4DynamicParticle(muon, G4ThreeVector()) );
  // ... finally, the parent
  products->SetParentParticle( *p );
  return products;
}


void G4DMuDFusionTChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "D_mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o << "The channel " << GetChannelName() 
      << " is only applicable to D_mu_D molecules, not "
      << part->GetParticleName();
    G4Exception("G4DMuDFusionTChannel::CheckIsApplicable",
                "DMUDF0003", FatalException, 
                o.str());

  }
}


// D_mu_D -> mu_T + n

G4DMuDFusionMuTChannel::
G4DMuDFusionMuTChannel(G4MuMolecule const* p, G4double fusionRate, 
			 G4int verboseLevel) :
  G4VMuMoleculeCaptureKineticsChannel(p, "D_mu_D -> mu_T + n fusion", 
				      fusionRate, verboseLevel),
  psdc(p->GetParticleName(), 1., 2, "neutron"){
  CheckIsApplicable();
  // mu_T
  psdc.SetDaughter(1,G4MuAtom::Definition(1,3));
}


G4DecayProducts* G4DMuDFusionMuTChannel::CaptureIt(G4DynamicParticle const* p){
  // FIXME ... not handling the molecular binding energy correctly
  // first, the fusion products ... The muon mass needs to be included
  // in the parent in the sticking case
  G4DecayProducts *products = psdc.DecayIt( p->GetMass() );
  // ... then, the parent
  products->SetParentParticle( *p );
  return products;
}


void G4DMuDFusionMuTChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "D_mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o << "The channel " << GetChannelName() 
      << " is only applicable to D_mu_D molecules, not "
      << part->GetParticleName();
    G4Exception("G4DMuDFusionMuTChannel::CheckIsApplicable",
                "DMUDF0004", FatalException, 
                o.str());
  }
}


// D_mu_D -> He4 + gamma + mu

G4DMuDFusionHe4Channel::G4DMuDFusionHe4Channel(G4MuMolecule const* p, G4double fusionRate, 
					       G4int verboseLevel) :
  G4VMuMoleculeCaptureKineticsChannel(p, "D_mu_D -> He4 + gamma + mu fusion", 
				      fusionRate, verboseLevel),
  psdc(p->GetParticleName(), 1., 2, "gamma", "alpha"){
  CheckIsApplicable();
}


G4DecayProducts* G4DMuDFusionHe4Channel::CaptureIt(G4DynamicParticle const* p){
  // FIXME ... not handling the molecular binding energy correctly
  G4MuonMinus *muon = G4MuonMinus::Definition();
  // first, the fusion products ...
  G4DecayProducts *products = 
    psdc.DecayIt( p->GetMass() - muon->GetPDGMass() );
  // ... then, the spectator muon ...
  products->PushProducts( new G4DynamicParticle(muon, G4ThreeVector()) );
  // ... finally, the parent
  products->SetParentParticle( *p );
  return products;
}


void G4DMuDFusionHe4Channel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "D_mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o << "The channel " << GetChannelName() 
      << " is only applicable to D_mu_D molecules, not "
      << part->GetParticleName();
    G4Exception("G4DMuDFusionHe4Channel::CheckIsApplicable",
                "DMUDF0005", FatalException, 
                o.str());    
  }
}


