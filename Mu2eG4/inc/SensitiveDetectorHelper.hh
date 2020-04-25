#ifndef Mu2eG4_SensitiveDetectorHelper_hh
#define Mu2eG4_SensitiveDetectorHelper_hh
//
// Some helper functions to manage repeated tasks related to sensitive detectors.
//
// $Id: SensitiveDetectorHelper.hh,v 1.5 2013/10/09 18:09:18 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/09 18:09:18 $
//
// Original author Rob Kutschke
//

// From Mu2e
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/Mu2eG4Config.hh"

// From the art tool chain
#include "cetlib/maybe_ref.h"

// From C++ and STL
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward references.
namespace art   { class Event; }
namespace art   { class ProducesCollector; }

namespace mu2e {

  class SimpleConfig;
  class SimParticleHelper;
  class Mu2eG4PerThreadStorage;

  class SensitiveDetectorHelper{

  public:

    SensitiveDetectorHelper(const Mu2eG4Config::SDConfig_& conf);
    // Accept compiler generator copy c'tor and assignment operator.

    // Register the sensitive detector with this class; to be called after G4 Initialize.
    void registerSensitiveDetectors();

    void declareProducts(art::ProducesCollector& collector);

    // Create data products and pre-fill with input hits if any;
    // to be called at the start of each event.
    void createProducts(const art::Event& evt, const SimParticleHelper& spHelper);

    // Attach the new per-event data products to the corresponding sensitive detector objects.
    void updateSensitiveDetectors(PhysicsProcessInfo& info,
                                  const SimParticleHelper& spHelper);

    // add the SD data into the PerThreadStorage
    void insertSDDataIntoPerThreadStorage(Mu2eG4PerThreadStorage* per_thread_store);

    //filter the event data here to cut down on execution time
    bool filterStepPointMomentum();
    bool filterTrackerStepPoints();

    // Query the same info
    bool enabled(StepInstanceName::enum_type instance) const;

    // Return one of the StepPointMCCollections.
    cet::maybe_ref<StepPointMCCollection> steps(StepInstanceName::enum_type id);

    // create SDs for arbitrary logical volumes as requested
    void instantiateLVSDs(const SimpleConfig& config);

    bool extMonPixelsEnabled() const { return extMonPixelsEnabled_; }
    ExtMonFNALPixelSD* getExtMonFNALPixelSD() const { return extMonFNALPixelSD_; }

    int verbosityLevel() const { return verbosityLevel_; }

  private:

    // A helper class to hold information about each sensitive detector object.
    struct StepInstance {
      explicit StepInstance(StepInstanceName::enum_type astepName):
        p(),
        stepName(StepInstanceName(astepName).name())
      {}

      explicit StepInstance(const std::string &name)  : stepName(name) {}

      explicit StepInstance():
        p(),
        stepName(){
      }

      // Accept compiler written d'tor, copy c'tor and assignment operator.

      // See note 2 in the .cc file for why the collection is not held by pointer.
      // For historical reasons the two names are different; maybe some day we will synchronize them.
      StepPointMCCollection    p;
      std::string              stepName;
      Mu2eSensitiveDetector *  sensitiveDetector = nullptr;
    };

    // Enabled pre-defined StepPointMC collections, except the timevd.
    typedef std::map<StepInstanceName::enum_type,StepInstance> InstanceMap;
    InstanceMap stepInstances_;

    // Logical volumes requested to be sensitive
    typedef std::map<std::string,StepInstance> LVSDMap;
    LVSDMap lvsd_;

    // Existing hit collections that should be merged into the output
    typedef std::vector<art::InputTag> InputTags;
    InputTags preSimulatedHits_;

    // Return all of the instances names of the data products to be produced.
    std::vector<std::string> stepInstanceNamesToBeProduced() const;

    // Separate handling as this detector does not produced StepPointMCs
    bool extMonPixelsEnabled_;
    ExtMonFNALPixelSD* extMonFNALPixelSD_ = nullptr;
    const bool standardMu2eDetector_;

    int  verbosityLevel_;

    // minimum momentum of a hit in a StepPtMCColl and minimum # hits in Tracker
    // to put event into art::Event
    double cutMomentumMin_;
    size_t minTrackerStepPoints_;
    std::vector<std::string> stepInstancesForMomentumCut_;

  };


} // end namespace mu2e
#endif /* Mu2eG4_SensitiveDetectorHelper_hh */
