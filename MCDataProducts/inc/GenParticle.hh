#ifndef MCDataProducts_GenParticle_hh
#define MCDataProducts_GenParticle_hh

//
// A temporary class to hold generated particles.
// It does not have a mother-daughter history.
//
//
// Original author Rob Kutschke
//
// I would prefer to use the HepMC/GenEvent and
// related classes but they are currently giving me
// compiler errors.
//
// In the mean time this class should be good enough.
//

// Mu2e includes
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenId.hh"

// Includes from external packages.
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

// C++ includes
#include <ostream>
#include <vector>

namespace mu2e {

  class GenParticle {

  public:

    // This c'tor is required for ROOT.
    GenParticle() : _time(), _properTime() {}

    GenParticle( PDGCode::type pdgId,
                 GenId generatorId,
                 CLHEP::Hep3Vector const&       position,
                 CLHEP::HepLorentzVector const& momentum,
                 double time,
                 double properTime = 0.):
      _pdgId(pdgId),
      _generatorId(generatorId),
      _position(position),
      _momentum(momentum),
      _time(time),
      _properTime(properTime)
    {}

    // Differs only in the type of second argument.
    GenParticle( PDGCode::type pdgId,
                 GenId::enum_type generatorId,
                 CLHEP::Hep3Vector const&       position,
                 CLHEP::HepLorentzVector const& momentum,
                 double time,
                 double properTime = 0.):
      _pdgId(pdgId),
      _generatorId(GenId(generatorId)),
      _position(position),
      _momentum(momentum),
      _time(time),
      _properTime(properTime)
    {}

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    // PDG particle ID code.
    PDGCode::type pdgId()       const { return _pdgId; }
    GenId         generatorId() const { return _generatorId; }
    double        time()        const { return _time;}

    // This is for multi-stage jobs, when we want to restart a
    // particle recorded in a previous job that already had non-zero
    // proper time.
    double        properTime()  const { return _properTime;}

    CLHEP::Hep3Vector const&       position() const { return _position;}
    CLHEP::HepLorentzVector const& momentum() const { return _momentum;}

    private:

    // PDG particle ID code.
    PDGCode::type _pdgId;

    // Identifier for which generator created this particle.
    GenId _generatorId;

    // Position and momentum at the given time.
    CLHEP::Hep3Vector       _position;
    CLHEP::HepLorentzVector _momentum;
    double _time;
    double _properTime;
  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const GenParticle& genp ){
    ost << "( "
        << genp.pdgId() << " "
        << genp.generatorId() << " "
        << genp.time() << " "
        << genp.properTime() << " "
        << genp.position() << " "
        << genp.momentum() << " "
        << " )";
    return ost;
  }
  typedef std::vector<mu2e::GenParticle> GenParticleCollection;
}

#endif /* MCDataProducts_GenParticle_hh */
