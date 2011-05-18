#ifndef ToyDP_GenParticle_hh
#define ToyDP_GenParticle_hh

//
// A temporary class to hold generated particles.
// It does not have a mother-daughter history.
//
// $Id: GenParticle.hh,v 1.1 2011/05/18 22:44:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 22:44:23 $
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
#include "ToyDP/inc/GenId.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// Includes from external packages.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {

  class GenParticle {

  public:

    // This c'tor is required for ROOT.
    GenParticle(){}

    GenParticle( PDGCode::type pdgId,
                    GenId generatorId,
                    CLHEP::Hep3Vector const&       position,
                    CLHEP::HepLorentzVector const& momentum,
                    double time):
      _pdgId(pdgId),
      _generatorId(generatorId),
      _position(position),
      _momentum(momentum),
      _time(time){
    }

    // Differs only in the type of second argument.
    GenParticle( PDGCode::type pdgId,
                 GenId::enum_type generatorId,
                 CLHEP::Hep3Vector const&       position,
                 CLHEP::HepLorentzVector const& momentum,
                 double time):
      _pdgId(pdgId),
      _generatorId(GenId(generatorId)),
      _position(position),
      _momentum(momentum),
      _time(time){
    }

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    // PDG particle ID code.
    PDGCode::type pdgId()       const { return _pdgId; }
    GenId         generatorId() const { return _generatorId; }
    double        time()        const { return _time;}

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

  };

}

#endif /* ToyDP_GenParticle_hh */
