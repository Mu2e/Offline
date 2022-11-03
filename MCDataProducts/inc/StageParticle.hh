// In some cases we want to re-sample an existing SimParticle and
// generate its daughters "by hand" rather than letting G4 do it.
// This class holds information produced by such custom generators
// and is used as input to Mu2eG4.
//
// Andrei Gaponenko, 2021

#ifndef MCDataProducts_StageParticle_hh
#define MCDataProducts_StageParticle_hh

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include <vector>
#include <ostream>

namespace mu2e {

  class StageParticle {
  public:

    StageParticle() = default;

    StageParticle(art::Ptr<SimParticle> parent,
                  ProcessCode creationCode,
                  PDGCode::type pdgId,
                  CLHEP::Hep3Vector position,
                  CLHEP::HepLorentzVector momentum,
                  double time)
      : parent_(parent)
      , creationCode_{creationCode}
      , pdgId_{pdgId}
      , position_{position}
      , momentum_{momentum}
      , time_{time}
    {}

    art::Ptr<SimParticle> parent() const { return parent_; }
    ProcessCode creationCode() const { return creationCode_; }
    PDGCode::type pdgId() const { return pdgId_; }
    double time() const { return time_;}
    const CLHEP::Hep3Vector& position() const { return position_;}
    const CLHEP::HepLorentzVector& momentum() const { return momentum_;}

  private:
    art::Ptr<SimParticle> parent_;
    ProcessCode creationCode_;
    PDGCode::type pdgId_ = PDGCode::unknown;
    CLHEP::Hep3Vector position_;
    CLHEP::HepLorentzVector momentum_;
    double time_ = 0;
  };

  typedef std::vector<mu2e::StageParticle> StageParticleCollection;

  std::ostream& operator<<(std::ostream& os, const StageParticle& s);
  std::ostream& operator<<(std::ostream& os, const StageParticleCollection& c);
}

#endif /* MCDataProducts_StageParticle_hh */
