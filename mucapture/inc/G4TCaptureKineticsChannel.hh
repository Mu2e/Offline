#ifndef G4TCaptureKineticsChannel_h
#define G4TCaptureKineticsChannel_h

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"

// FIXME ... should I add an IsApplicable member and a check on
// ckt::Insert?  That would be nice of me :-)

template<class T> class G4TCaptureKineticsChannel {

public:

  G4TCaptureKineticsChannel(T const* particle,
			    G4String const& channelName, 
			    G4double captureRate,
			    G4int verboseLevel = 0);

  virtual ~G4TCaptureKineticsChannel(){};

  virtual G4bool IsCloneable() const { return false; }
  virtual G4TCaptureKineticsChannel* Clone(T const*){ return 0; }

  G4double GetRate() const { return rate; }

  virtual G4DecayProducts* CaptureIt(G4DynamicParticle const*) = 0;

  void SetVerboseLevel(G4int level) { verbose = level; }
  G4int GetVerboseLevel() const { return verbose; }

  G4String GetChannelName() const { return name; }

protected:
  G4TCaptureKineticsChannel(G4TCaptureKineticsChannel const&);
  G4TCaptureKineticsChannel& operator=(G4TCaptureKineticsChannel const&);

  T const* part;
  G4String name;
  G4double rate;
  G4int verbose;

private:

};

#define G4TCAPTUREKINETICSCHANNEL_ICC
#include "G4TCaptureKineticsChannel.icc"
#undef G4TCAPTUREKINETICSCHANNEL_ICC

#endif
