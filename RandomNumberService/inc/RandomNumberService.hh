#ifndef RandomNumberService_RandomNumberService_hh
#define RandomNumberService_RandomNumberService_hh
//
// Maintain multiple independent random numbers engines.
// This includes saving and restoring seeds and state.
//
// $Id: RandomNumberService.hh,v 1.2 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/05 16:07:38 $
//
// Original author Rob Kutschke
//
// Notes:
//  1) The present implementation just seeds the CLHEP global
//     random number generator, HepRandom.
//
//  2) For details on the planned future behaviour:
//       http://mu2e.fnal.gov/atwork/computing/Random.shtml
//
//  3) These methods are examples of the available callbacks.
//     I was experimenting to check which ones are useful.
//

// C++ include files
#include <string>

// Framework include files
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

// CLHEP include files.
#include "CLHEP/Random/RandomEngine.h"

// Forward declarations.
namespace edm{
  class ParameterSet;
  class ActivityRegistry;
  class EventId;
  class Event;
  class EventSetup;
  class Timestamp;
}

namespace mu2e {

  class RandomNumberService : public edm::RandomNumberGenerator {

  public:
    typedef std::vector<std::string>            LabelInfo;
    typedef std::vector<std::vector<uint32_t> > StateInfo;
    typedef std::vector<std::vector<uint32_t> > SeedInfo;
    
    RandomNumberService(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~RandomNumberService();
    
    // See note 3.
    void postBeginJob();
    void postEndJob();
    void preSource();
    void postSource();
    void preProcessEvent(edm::EventID const& id, edm::Timestamp const& iTime);
    void postProcessEvent(edm::Event const&, edm::EventSetup const&);
 
    // The methods below are required by the RandomNumberGenerator interface
    
    // This  methods is the main method called by user code.
    // Also called by IOPool/Input/src/RootInputFileSequence.cc
    CLHEP::HepRandomEngine& getEngine() const;

    // Obsolete.  Should never be called and will go away.
    uint32_t mySeed() const;
    
    // Called by RandomNumberService/src/RandomNumberSaver_plugin.cc
    const LabelInfo& getCachedLabels() const {return _labels; }
    const StateInfo& getCachedStates() const {return _states; }
    const SeedInfo&  getCachedSeeds()  const {return _seeds;  }

    // Called by FWCore/Framework/src/InputSource.cc
    void snapShot();
    void restoreState(const edm::Event& event);

    // Will be called by RandomNumberSaver_plugin.cc or maybe by
    // by one of the callback methods?
    void saveEngineState(const std::string& fileName);
    void restoreEngineState(const std::string& fileName);

    // For debugging purposes only
    void print();

    // end of methods required by RandomNumberGenerator

  private:
    
    // Temporary hack.  See note 1.
    long _globalSeed;

    // The per-module-instance information.
    LabelInfo _labels;
    StateInfo _states;
    SeedInfo  _seeds;

  };
}

#endif
