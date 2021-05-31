// Reject events with incorrect hits (initially, outside the given detector volume)
// Heavily based on ReadBack_module
// Krzysztof Genser, 2015

#include <string>
#include <iostream>
#include <cmath>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e {

  //================================================================
  class VetoIncorrectHits : public art::EDFilter {
  public:
    explicit VetoIncorrectHits(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
    bool doTracker(const art::Event& event);
    bool doCRV(const art::Event& event);

  private:

    // Diagnostics printout level
    int diagLevel_;

    // Limit on number of events for which there will be full printout.
    int maxFullPrint_;

    // sould throw or not
    bool doNotThrow_;

    // Module label of the geant4 module that made the hits.
    std::string g4ModuleLabel_;

    // Name of the tracker StepPoint collection
    std::string trackerStepPoints_;


    // by how much the straw gas and hit relative positions can differ
    double strawHitPositionTolerance_;

    // Name of the CRSScintillatorBar(CRV) StepPoint collection
    std::string crvStepPoints_;

    // by how much the crv bar and hit relative positions can differ
    double crvHitPositionTolerance_;

    // End: run time parameters

    // Number of events processed.
    int nProcessed_;

    // a map of passed events per object
    std::map<std::string,int> passedEventsMap_;

    // a map of seen hits per object
    std::map<std::string,int> seenHitsMap_;

    // a map of passed hits per object
    std::map<std::string,int> passedHitsMap_;


  };

  //================================================================
  VetoIncorrectHits::VetoIncorrectHits(fhicl::ParameterSet const& pset)
    :
    art::EDFilter{pset},
    // Run time parameters
    diagLevel_(pset.get<int>("diagLevel")),
    maxFullPrint_(pset.get<int>("maxFullPrint")),
    doNotThrow_(pset.get<bool>("doNotThrow")),
    g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
    trackerStepPoints_(pset.get<std::string>("trackerStepPoints")),
    strawHitPositionTolerance_(pset.get<double>("shPositionTolerance")),
    crvStepPoints_(pset.get<std::string>("crvStepPoints")),
    crvHitPositionTolerance_(pset.get<double>("crvPositionTolerance")),
    // End of run time parameters
    nProcessed_(0),
    passedEventsMap_(std::map<std::string,int>()),
    passedHitsMap_(std::map<std::string,int>())
  {
    std::cout << __func__ << std::endl << pset.to_indented_string();
  }
  //================================================================
  bool VetoIncorrectHits::filter(art::Event& theEvent) {

    // will not modify the event, so force the const
    const art::Event& event = const_cast<art::Event&>(theEvent);
    bool passed = true;

    bool trackerResult =  doTracker(event);
    bool crvResult     =  doCRV(event);

    // a simple && for now

    passed = trackerResult && crvResult;

    ++nProcessed_;

    return passed;

  }

  //================================================================
  bool VetoIncorrectHits::doCRV(const art::Event& event){

    bool passed = true;

    // Get a reference to CosmicRayShield (it contains crv)
    GeomHandle<CosmicRayShield> cosmicRayShieldGeomHandle;
    CosmicRayShield const& crv(*cosmicRayShieldGeomHandle);

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(g4ModuleLabel_,crvStepPoints_,hits);
    // no hits means hits are "good"
    if ( ! hits.isValid() ) { return passed; }

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      ++seenHitsMap_[crvStepPoints_];

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();

      // Get the CRSScintillatorBar information:
      const CRSScintillatorBar&  bar = crv.getBar( hit.barIndex() );
      CLHEP::Hep3Vector const &  mid = bar.getPosition();

      CLHEP::Hep3Vector hitLocal  = pos-mid;
      CLHEP::Hep3Vector hitLocalN = CLHEP::Hep3Vector(hitLocal.x()/bar.getHalfLengths()[0],
                                                      hitLocal.y()/bar.getHalfLengths()[1],
                                                      hitLocal.z()/bar.getHalfLengths()[2]);

      // Print out limited to the first few events.
      if (  diagLevel_ > 1 && nProcessed_ < maxFullPrint_ ){

        // mf::LogInfo("GEOM") << "VetoIncorrectHits::doCRV"
        std::cout << "VetoIncorrectHits::doCRV"
                  << " Event/Hit: "
                  << event.id().event() << " as "
                  << nProcessed_ << " "
                  << i                  <<  " "
                  << hit.trackId()      << "   "
                  << hit.volumeId()     << " "
                  << bar.id()           << " | "
                  << pos                << " "
                  << mid                << " "
                  << (mid-pos)          << " | "
                  << hitLocalN          << " | "
                  << std::endl;

      }

      if ( ( std::abs(hitLocalN.x()) - 1. > crvHitPositionTolerance_ ) ||
           ( std::abs(hitLocalN.y()) - 1. > crvHitPositionTolerance_ ) ||
           ( std::abs(hitLocalN.z()) - 1. > crvHitPositionTolerance_ )
           ) {

        passed = false;

        std::ostringstream os;
        os << __func__
           << " Event " << event.id().event() << " as " << nProcessed_
           << " Hit " << i << " "
           << pos <<  " "
           << " ouside the crv bar " << bar.id()
           << " inconsistent tracker geometry file? "
           << "; relative differences: "
           << std::abs(hitLocalN.x()) - 1. << ", "
           << std::abs(hitLocalN.y()) - 1. << ", "
           << std::abs(hitLocalN.z()) - 1.
           << "; tolerance : "
           << crvHitPositionTolerance_
           << std::endl;

        if ( doNotThrow_ ) {

          if ( diagLevel_ > -1 &&  nProcessed_< maxFullPrint_ ){
            // mf::LogWarning("GEOM") << os.str();
            std::cout << os.str();
          }

        } else {
          throw cet::exception("GEOM") << os.str();
        }

      }

      passed && ++passedHitsMap_[crvStepPoints_];

    } // end loop over hits.

    passed && ++passedEventsMap_[crvStepPoints_];
    return passed;

  } // end doCRV

  //================================================================
  bool VetoIncorrectHits::doTracker(const art::Event& event){

    // Get a reference to the tracker
    // Throw exception if not successful.
    const Tracker& tracker = *GeomHandle<Tracker>();

    bool passed = true;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(g4ModuleLabel_,trackerStepPoints_,hits);

    // no hits means hits are "good"
    if ( ! hits.isValid() ) { return passed; }

    // Loop over all hits.
    for ( size_t i=0; i<hits->size(); ++i ){

      ++seenHitsMap_[trackerStepPoints_];

      // Alias, used for readability.
      const StepPointMC& hit = (*hits)[i];

      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get the straw information:
      const Straw&   straw = tracker.getStraw( hit.strawId() );
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);

      // Get the radius of the reference point in the local
      // coordinates of the straw
      double s = w.dot(pos-mid);
      CLHEP::Hep3Vector point = pos - (mid + s*w);

      double normPointMag = point.mag()/tracker.strawInnerRadius();
      double normS = s/straw.halfLength();

      if ( diagLevel_ > 1 &&  nProcessed_< maxFullPrint_ ){

        std::cout << "VetoIncorrectHits::doTracker"
                  << " Event " << event.id().event() << " as " << nProcessed_
                  << " normalized reference point - 1 : "
                  << std::scientific
                  << normPointMag - 1.
                  << " normalized wire z of reference point - 1 : "
                  << std::abs(normS) -1.
                  << std::fixed
                  << std::endl;

      }

      if ( ( normPointMag - 1. > strawHitPositionTolerance_ ) ||
           ( std::abs(normS) - 1. > strawHitPositionTolerance_ )
           ) {

        passed = false;

        std::ostringstream os;
        os << __func__
           << " Event " << event.id().event() << " as " << nProcessed_
           << " Hit " << i << " "
           << pos <<  " "
           << " ouside the straw " << straw.id()
           << " inconsistent tracker geometry file? "
           << "; radial difference: "
           << normPointMag - 1.
           << ", longitudinal difference: "
           << std::abs(normS) - 1.
           << "; tolerance : "
           << strawHitPositionTolerance_
           << std::endl;

        if ( doNotThrow_ ) {

          if ( diagLevel_ > -1 &&  nProcessed_< maxFullPrint_ ){
            // mf::LogWarning("GEOM") << os.str();
             std::cout << os.str();
          }

        } else {
          throw cet::exception("GEOM") << os.str();
        }

      }

      passed && ++passedHitsMap_[trackerStepPoints_];

    } // end loop over hits.

    passed && ++passedEventsMap_[trackerStepPoints_];

    return passed;

  }

  //================================================================
  void VetoIncorrectHits::endJob() {

    std::ostringstream os;
    os << "VetoIncorrectHits_module "
       << std::setw(20)
       << "summary:  processed  :"
       << " "
       << std::setw(20)
       << nProcessed_ << " events" << std::endl;

    //    std::cout << " VetoIncorrectHits::endJob size of the map " << statMap_.size() << std::endl;

    for(const auto &i : passedEventsMap_) {
      os << "VetoIncorrectHits_module "
         << std::setw(20) << i.first << " : "
         << std::setw(20) << i.second
         << " events passed  " << std::endl
         << "VetoIncorrectHits_module "
         << std::setw(20) << i.first << " : "
         << std::setw(20) << seenHitsMap_[i.first]
         << " hits processed " << std::endl
         << "VetoIncorrectHits_module "
          << std::setw(20) << i.first << " : "
         << std::setw(20) << passedHitsMap_[i.first]
         << " hits passed    " << std::endl;
    }

    mf::LogInfo("Summary")<<os.str();

  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VetoIncorrectHits);
