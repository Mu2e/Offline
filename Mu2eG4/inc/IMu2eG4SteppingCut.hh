// The interface for Mu2eG4SteppingAction cuts.
//
// Andrei Gaponenko, 2015

#ifndef Mu2eG4_IMu2eG4SteppingCut_hh
#define Mu2eG4_IMu2eG4SteppingCut_hh

class G4Step;
namespace CLHEP { class Hep3Vector; }

namespace art { class Event; }
namespace art { class EDProducer; }

namespace mu2e {

  class SimParticleHelper;

  class IMu2eG4SteppingCut {
  public:

    virtual bool evaluate(const G4Step *aStep) = 0;

    virtual void declareProducts(art::EDProducer *parent) =  0;

    virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) = 0;

    // Create data  products and pre-fill with input hits, if any; to be called at the start of each event.
    virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) = 0;

    // Put the data products into the event.
    virtual void put(art::Event& event) = 0;

    virtual ~IMu2eG4SteppingCut() {}
  };

} // end namespace mu2e

#endif /* Mu2eG4_IMu2eG4SteppingCut_hh */
