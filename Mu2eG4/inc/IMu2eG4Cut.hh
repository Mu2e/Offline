// This file defines an interface for making cuts in
// Mu2eG4SteppingAction and Mu2eG4StackingAction.  Example cuts are a
// cut on the G4 volume of a track/step (aka
// "g4SteppingAction.killInTheseVolumes" in G4_module's SimpleConfig),
// or a cut on the energy of a particle ("g4.killLowEKine" in
// SimpleConfig).  Unlike the SimpleConfig situation we can combine
// Mu2eG4_module cuts in a flexible way using fcl.  For example, we
// can kill low energy particles only in specific volumes.  New cuts
// include a cut on the position of a particle in the Mu2e coordinate
// system, independently of G4 volumes.

//
// In addition to providing a cut decision to the stepping action, the
// code that evaluates a stepping cut can write information into a
// data product.  For multi-stage job setups we will be stopping
// tracks and writing out StepPointMCs, which will be used to re-start
// the tracks in a subsequent stage.  Note that geometric cuts, which
// are independent of G4 volumes, will triger in the steping action
// after the track has propagated past the cut surface by some amount.
// Because a purpose of a multi-stage cut is to hide details of the
// apparatus behind the cut surface from the that stage of
// simulations, we will write out the *beginning* of the step when a
// cut condition is satisfied at the *end* point of a step.
//
// Cuts in the stacking action SHOULD NOT write out info.  At least
// not into the same collection as the stepping cuts.  Otherwise we
// can end up with double counting: a step which is being cut and
// saved in the stepping action can produce a secondary at its end,
// which would be collected by the stacking cut and re-used by the
// next stage despite of the last steb being re-simulated as well.
//
// Andrei Gaponenko, 2015

#ifndef Mu2eG4_IMu2eG4Cut_hh
#define Mu2eG4_IMu2eG4Cut_hh

#include <memory>

class G4Step;
class G4Track;

namespace CLHEP { class Hep3Vector; }

namespace art { class Event; }
namespace art { class ProducesCollector; }
namespace fhicl { class ParameterSet; }

namespace mu2e {

  class SimParticleHelper;
  class Mu2eG4ResourceLimits;
  class Mu2eG4PerThreadStorage;

  class IMu2eG4Cut {
  public:

    virtual bool steppingActionCut(const G4Step  *step) = 0;
    virtual bool stackingActionCut(const G4Track *trk) = 0;

    virtual void declareProducts(art::ProducesCollector& collector) =  0;

    virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) = 0;

    // Create data products and pre-fill with input hits, if any; to be called at the start of each event.
    virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) = 0;

    // put data into the stash
    virtual void insertCutsDataIntoPerThreadStorage(Mu2eG4PerThreadStorage* per_thread_store) = 0;

    // delete data if we don't need it
    virtual void deleteCutsData() = 0;

    // Put the data products into the event.
    //virtual void put(art::Event& event) = 0;

    virtual ~IMu2eG4Cut() {}
  };

  //================================================================
  std::unique_ptr<IMu2eG4Cut> createMu2eG4Cuts(const fhicl::ParameterSet& pset, const Mu2eG4ResourceLimits& mu2elimits);

} // end namespace mu2e

#endif /* Mu2eG4_IMu2eG4Cut_hh */
