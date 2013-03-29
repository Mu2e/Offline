#ifndef Mu2eG4_SensitiveDetectorHelper_hh
#define Mu2eG4_SensitiveDetectorHelper_hh
//
// Some helper functions to manage repeated tasks related to sensitive detectors.
//
// $Id: SensitiveDetectorHelper.hh,v 1.2 2013/03/29 04:35:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/03/29 04:35:17 $
//
// Original author Rob Kutschke
//

// From Mu2e
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

// From the art tool chain
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/maybe_ref.h"

// From C++ and STL
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward references.
namespace art   { class Event; }
namespace fhicl { class ParameterSet; }

namespace mu2e {

  class SensitiveDetectorHelper{

  public:

    SensitiveDetectorHelper( fhicl::ParameterSet const& pset );
    // Accept compiler generator d'tor, copy c'tor and assignment operator.

    // Register the sensitive detector with this class; to be called after G4 Initialize.
    void registerSensitiveDetectors();

    // Create empty data products; to be called at the start of each event.
    void createProducts();

    // Attach the new per-event data products to the corresponding sensitive detector objects.
    void updateSensitiveDetectors( PhysicsProcessInfo&   info,
                                   art::ProductID const& simsId,
                                   art::Event&           event );

    // Put the data products into the event.
    void put( art::Event& event);

    // Return all of the instances names of the data products to be produced.
    std::vector<std::string> stepInstanceNamesToBeProduced() const;

    // Query the same info
    bool enabled(StepInstanceName::enum_type instance) const;

    // Return one of the StepPointMCCollections.
    cet::maybe_ref<StepPointMCCollection> steps( StepInstanceName::enum_type id );

  private:

    // A helper class to hold information about each sensitive detector object.
    struct StepInstance {
      explicit StepInstance(StepInstanceName::enum_type astepName ):
        p(),
        stepName(StepInstanceName(astepName)),
        sensitiveDetector(0){
      }

      explicit StepInstance():
        p(),
        stepName(),
        sensitiveDetector(0){
      }

      // Accept compiler written d'tor, copy c'tor and assignment operator.

      // See note 1 in the .cc file for why the collection is not held by pointer.
      // For historical reasons the two names are different; maybe some day we will synchronize them.
      StepPointMCCollection    p;
      StepInstanceName         stepName;
      Mu2eSensitiveDetector *  sensitiveDetector;

    };

    // All of the StepPointMC collections, except the timevd.
    //std::vector<StepInstance> stepInstances_;
    typedef std::map<StepInstanceName::enum_type,StepInstance> InstanceMap;
    InstanceMap stepInstances_;

  };


} // end namespace mu2e
#endif /* Mu2eG4_SensitiveDetectorHelper_hh */
