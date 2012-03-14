#ifndef G4DMuDFusionChannels_hh
#define G4DMuDFusionChannels_hh 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 30 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuMolecule.hh"
#include "G4VMuMoleculeCaptureKineticsChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

// There are five (or perhaps six?) fusion channels:

// 1) DmuD(+) -> He3(++) + n + mu-
// 2)         -> muHe3(+) + n
// 3)         -> T(+) + T(+) + mu(-)
// 4)         -> muT + p(+)
// 5)         -> He4(++) + gamma + mu(-)
// 6) ???     -> muHe4(+) + gamma

// Prerequisites:  Be sure to instantiate the following in your
// physics list:
// G4MuonMinus, G4Gamma, G4Proton, G4Neutron, G4He3, G4Alpha,
// G4Triton

class G4DMuDFusionHe3Channel : public G4VMuMoleculeCaptureKineticsChannel {
public:
  G4DMuDFusionHe3Channel(G4MuMolecule const* p, G4double fusionRate, 
			 G4int verboseLevel = 0);
  G4DecayProducts *CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4DMuDFusionMuHe3Channel : public G4VMuMoleculeCaptureKineticsChannel {
public:
  G4DMuDFusionMuHe3Channel(G4MuMolecule const* p, G4double fusionRate, 
			   G4int verboseLevel = 0);
  G4DecayProducts *CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4DMuDFusionTChannel : public G4VMuMoleculeCaptureKineticsChannel {
public:
  G4DMuDFusionTChannel(G4MuMolecule const* p, G4double fusionRate, 
		       G4int verboseLevel = 0);
  G4DecayProducts *CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4DMuDFusionMuTChannel : public G4VMuMoleculeCaptureKineticsChannel {
public:
  G4DMuDFusionMuTChannel(G4MuMolecule const* p, G4double fusionRate, 
			 G4int verboseLevel = 0);
  G4DecayProducts *CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};

class G4DMuDFusionHe4Channel : public G4VMuMoleculeCaptureKineticsChannel {
public:
  G4DMuDFusionHe4Channel(G4MuMolecule const* p, G4double fusionRate, 
			 G4int verboseLevel = 0);
  G4DecayProducts *CaptureIt(G4DynamicParticle const*);
private:
  void CheckIsApplicable() const;
  G4PhaseSpaceDecayChannel psdc;
};



#endif
