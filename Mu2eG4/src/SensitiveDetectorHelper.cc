//
// A helper class manage repeated tasks related to sensitive detectors.
//
//
// Original author Rob Kutschke
//
//
// Notes.
//
// 1) The timeVD virtual detector does not use the G4 sensitve detector methodology.
//    It is, instead, implemented inside SteppingAction.  So it is not managed
//    by this class.
//
// 2) The embedded class StepIntance needs to be copyable so that it can be put into
//    the std::vector.  Therefore it cannot contain an unique_ptr to its StepPointMC
//    collection.  The work around is to hold the collection by value and to use swap
//    to transfer it into the unique_ptr that will be given to the event.  This is
//    a very small CPU time penalty but it saves us from doing any explicit memory management.
//

// From Mu2e
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "GeometryService/inc/GeometryService.hh"

// From art and its tool chain
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// From G4
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Threading.hh"

#include <map>

using namespace std;

namespace mu2e {

  //================================================================
  SensitiveDetectorHelper::SensitiveDetectorHelper(const Mu2eG4Config::SDConfig_& conf)
    :
    extMonPixelsEnabled_(false),
    standardMu2eDetector_((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    verbosityLevel_(conf.verbosityLevel()),
    cutMomentumMin_(conf.cutMomentumMin()),
    minTrackerStepPoints_(conf.minTrackerStepPoints())
  {

    bool enableAllSDs = false;
    conf.enableAllSDs(enableAllSDs);

    std::vector<std::string> enabledSDs;
    if(conf.enableSD(enabledSDs)) {
      if(enableAllSDs) {
        throw cet::exception("CONFIG")<<"SensitiveDetectorHelper: enabledSDs should not be used when enableAllSDs=true\n";
      }
    }
    else {
      if(!enableAllSDs) {
        throw cet::exception("CONFIG")<<"SensitiveDetectorHelper: enabledSDs must be specified used when enableAllSDs is not set to true\n";
      }
    }

    typedef std::set<std::string> TODO;
    TODO todo(enabledSDs.begin(), enabledSDs.end());

    if(enableAllSDs || (todo.find(SensitiveDetectorName::ExtMonFNAL()) != todo.end())) {
      extMonPixelsEnabled_ = true;
      todo.erase(SensitiveDetectorName::ExtMonFNAL());
    }

    // Build list of StepInstances to look after.  See note 1.
    std::vector<StepInstanceName> const& preDefinedSD(StepInstanceName::allValues());
    for ( std::vector<StepInstanceName>::const_iterator i=preDefinedSD.begin(); i != preDefinedSD.end(); ++i ){

      if(enableAllSDs || (todo.find(i->name()) != todo.end())) {
        todo.erase(i->name());

        // stepper does not have a sensitive detector associated with it...
        if (*i == StepInstanceName::stepper) continue;
        StepInstanceName::enum_type id(*i);
        if ( id != StepInstanceName::unknown && id !=StepInstanceName::timeVD ){
          stepInstances_[id] = StepInstance(id);
        }//if
      }//if enableAllSDs
    }//for

    if(!todo.empty()) {
      std::ostringstream os;
      os << "SensitiveDetectorHelper: the following detectors in the enableSD list are not recognized: ";
      std::copy(todo.begin(), todo.end(), std::ostream_iterator<std::string>(os, " "));
      throw cet::exception("CONFIG")<<os.str();
    }

    //----------------
    std::vector<string> lvlist(conf.sensitiveVolumes());
    for(const auto& name : lvlist) {
      lvsd_[name] = StepInstance(name);
    }

    //----------------
    // Careful here: we are in the constructor, so class methods may not be used freely.
    // However calling the following at OK at this point (but not earlier!)

    const vector<string> outputs(stepInstanceNamesToBeProduced());

    //----------------
    // New hits will be added to hits from these existing collections
    for(const auto& s : conf.preSimulatedHits()) {
      preSimulatedHits_.emplace_back(s);
      // check that the destination is known and enabled
      if(std::find(outputs.begin(), outputs.end(), preSimulatedHits_.back().instance()) == outputs.end()) {
        throw cet::exception("CONFIG")<<"SensitiveDetectorHelper: no matching destination for preSimulatedHits = "<<s<<"\n";
      }//if
    }//for

    //----------------

    for(const auto& i : conf.inputs()) {
      stepInstancesForMomentumCut_.emplace_back(i);
    }

  }//end c'tor

  //================================================================
  void SensitiveDetectorHelper::instantiateLVSDs(const SimpleConfig& config){

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    art::ServiceHandle<Mu2eG4Helper> helper;

    for(auto& iter : lvsd_) {
      iter.second.sensitiveDetector = new Mu2eG4SensitiveDetector(iter.first, config);
      SDman->AddNewDetector(iter.second.sensitiveDetector);
      helper->locateVolInfo(iter.first).logical->SetSensitiveDetector(iter.second.sensitiveDetector);
    }
  }


  // Find the sensitive detector objects and attach them to each StepInstance object.
  // Must not be called until G4 has been initialized.
  void SensitiveDetectorHelper::registerSensitiveDetectors(){
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ) {
      StepInstance& step(i->second);
            if (verbosityLevel_ > 0 ) {
        std::cout << __func__ << " looking for sd with name: "
                  << step.stepName << std::endl;
            }
      bool printWarnings = (verbosityLevel_ > -1) ? true : false;
      step.sensitiveDetector =
        dynamic_cast<Mu2eG4SensitiveDetector*>(sdManager->FindSensitiveDetector(step.stepName.c_str(),printWarnings));
    }

    extMonFNALPixelSD_ = ( standardMu2eDetector_ && extMonPixelsEnabled_) ?
      dynamic_cast<ExtMonFNALPixelSD*>(sdManager->
                                       FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()))
      : nullptr;

  }


  // Create new data products.  To be called at start of each event.
  void SensitiveDetectorHelper::createProducts(const art::Event& event,
                                               const SimParticleHelper& spHelper){
    //----------------
    // Read in pre-simulated hits that we want to merge into the outputs

    // instance => multiple hit collections
    typedef std::multimap<string, art::ValidHandle<StepPointMCCollection> > InHits;
    InHits inputHits;

    for(const auto& tag : preSimulatedHits_) {
      inputHits.insert(make_pair(tag.instance(), event.getValidHandle<StepPointMCCollection>(tag)));
    }

    //----------------
    // Clean and pre-fill pre-defined SD collections in stepInstances_

    for (auto& i : stepInstances_) {

      auto& out = i.second.p;
      out.clear();

      // Copy all input collection with the current instance name
      const auto rr = inputHits.equal_range(i.second.stepName);
      for(auto in = rr.first; in != rr.second; ++in) {
        out.insert(out.end(), in->second->cbegin(), in->second->cend());
      }

      // Update SimParticle ptr to point to the new collection
      for(auto& hit : out) {
        hit.simParticle() =
          art::Ptr<SimParticle>(spHelper.productID(),
                                hit.simParticle()->id().asUint(),
                                spHelper.productGetter());

      }//for auto& hit
    }//for auto& i

    //----------------
    // Clean and pre-fill the logical volume collections

    for(auto& i : lvsd_) {
      auto& out = i.second.p;
      out.clear();

      // Copy all input collection with the current instance name
      const auto rr = inputHits.equal_range(i.first);
      for(auto in = rr.first; in != rr.second; ++in) {
        out.insert(out.end(), in->second->cbegin(), in->second->cend());
      }

      // Update SimParticle ptr to point to the new collection
      for(auto& hit : out) {
        hit.simParticle() =
          art::Ptr<SimParticle>(spHelper.productID(),
                                hit.simParticle()->id().asUint(),
                                spHelper.productGetter());
      }//for auto& hit
    }//for auto& i
  }


  // Hand each sensitive detector object the StepPointMC collection that it will fill.
  // Also hand it the information needed to create art::Ptr objects.
  void SensitiveDetectorHelper::updateSensitiveDetectors( PhysicsProcessInfo&   info,
                                                          const SimParticleHelper& spHelper){

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ){

      StepInstance& instance(i->second);
      if ( instance.sensitiveDetector ){
        instance.sensitiveDetector->beforeG4Event(instance.p, info, spHelper);
      }//if
    }//for

    for(auto& i : lvsd_) {
      i.second.sensitiveDetector->beforeG4Event(i.second.p, info, spHelper);
    }//for

  }


  void SensitiveDetectorHelper::insertSDDataIntoPerThreadStorage(Mu2eG4PerThreadStorage* per_thread_store){

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ) {
      unique_ptr<StepPointMCCollection> p(new StepPointMCCollection);
      StepInstance& instance(i->second);
      std::swap( instance.p, *p);
      per_thread_store->insertSDStepPointMC(std::move(p), instance.stepName);
    }

    for (auto& i: lvsd_) {
      unique_ptr<StepPointMCCollection> p(new StepPointMCCollection);
      std::swap( i.second.p, *p);
      per_thread_store->insertSDStepPointMC(std::move(p), i.second.stepName);
    }
  }


  bool SensitiveDetectorHelper::filterStepPointMomentum(){

    bool passed = false;

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ) {
      StepInstance& instance(i->second);

      for (std::vector<std::string>::iterator j=stepInstancesForMomentumCut_.begin();
           j != stepInstancesForMomentumCut_.end(); j++) {

        if (instance.stepName == *j) {
          for(const auto& hit : instance.p) {
            if(hit.momentum().mag() > cutMomentumMin_) {
              passed = true;
              break;
            }//if momentum
          }//for hit
        }//if stepName
      }//for stepInstances
    }

    for (auto& i: lvsd_) {

      for (std::vector<std::string>::iterator j=stepInstancesForMomentumCut_.begin();
           j != stepInstancesForMomentumCut_.end(); j++) {

        if (i.second.stepName == *j) {

          for(const auto& hit : i.second.p) {
            if(hit.momentum().mag() > cutMomentumMin_) {
              passed = true;
              break;
            }//if momentum
          }//for hit
        }//if stepName
      }//for stepInstancesForMomentumCut_
    }

    //++numInputEvents_;
    //if(passed) { ++numPassedEvents_; }
    //void FilterStepPointMomentum::endJob() {
    //    mf::LogInfo("Summary")
    //    <<"FilterStepPointMomentum_module: passed "
    //    <<numPassedEvents_<<" / "<<numInputEvents_<<" events\n";
    //}

    return passed;
  }


  bool SensitiveDetectorHelper::filterTrackerStepPoints(){

    bool passed = false;

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ) {
      if (i->second.stepName == "tracker" && i->second.p.size() >= minTrackerStepPoints_) {
        passed = true;
      }
    }//for stepInstances

    for (auto& i: lvsd_) {

      for (std::vector<std::string>::iterator j=stepInstancesForMomentumCut_.begin();
           j != stepInstancesForMomentumCut_.end(); j++) {

        if (i.second.stepName == "tracker" && i.second.p.size() >= minTrackerStepPoints_) {
          passed = true;
        }

      }
    }
    return passed;
  }

  // Return all of the instances names of the data products to be produced.
  vector<string> SensitiveDetectorHelper::stepInstanceNamesToBeProduced() const{

    vector<string> names;
    for ( InstanceMap::const_iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ){
      names.push_back( i->second.stepName );
    }
    for(const auto& i : lvsd_) {
      names.emplace_back(i.second.stepName);
    }

    return names;
  }


  void SensitiveDetectorHelper::declareProducts(art::ProducesCollector& collector) {

    vector<string> const& instanceNames = stepInstanceNamesToBeProduced();
    for(const auto& name: instanceNames) {
      collector.produces<StepPointMCCollection>(name);
    }
    if(extMonPixelsEnabled_)
      collector.produces<ExtMonFNALSimHitCollection>();
  }


  bool SensitiveDetectorHelper::enabled(StepInstanceName::enum_type instance) const{
    return stepInstances_.find(instance) != stepInstances_.end();
  }


  // Return one of the StepPointMCCollections.
  cet::maybe_ref<StepPointMCCollection> SensitiveDetectorHelper::steps( StepInstanceName::enum_type id ){

    InstanceMap::iterator i = stepInstances_.find(id);
    if ( i != stepInstances_.end() ) {
      return cet::maybe_ref<StepPointMCCollection>(i->second.p);
    }
    // Cannot find the requested id; return an empty maybe_ref.
    return cet::maybe_ref<StepPointMCCollection>();
  }

}  // end namespace mu2e
