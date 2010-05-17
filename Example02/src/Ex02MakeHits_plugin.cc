/*----------------------------------------------------------------------

  A producer module that makes a collection of overly simplified "hits"
  and adds them to the event.

  $Id: Ex02MakeHits_plugin.cc,v 1.4 2010/05/17 21:47:33 genser Exp $
  $Author: genser $
  $Date: 2010/05/17 21:47:33 $
   
  Original author Rob Kutschke


----------------------------------------------------------------------*/

// C++ includes.
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Other external includes.
#include <boost/shared_ptr.hpp>

// Mu2e includes.
#include "ToyDP/inc/ToyHitCollection.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CTrackerGeom/inc/CTracker.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "GeneralUtilities/inc/pow.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class Ex02MakeHits : public edm::EDProducer {
    
  public:
    explicit Ex02MakeHits(edm::ParameterSet const& pSet) : 
      minPulseHeight_(pSet.getParameter<double>("minPulseHeight"))
    {
      produces<ToyHitCollection>();
    }
    virtual ~Ex02MakeHits() { }

    virtual void produce(edm::Event& e, edm::EventSetup const& c);
    
  private:

    // Parameters from run time configuration.
    double minPulseHeight_;

  };

  void
  Ex02MakeHits::produce(edm::Event& event, edm::EventSetup const&) {

    // Geometry information about the CTracker.
    GeomHandle<CTracker> ctracker;

    // Extract the event number from the event.
    int eventNumber = event.id().event();

    // Make the collection to hold the hits that we will make.
    auto_ptr<ToyHitCollection> p(new ToyHitCollection);

    // Instantiate an object that can generate random vectors uniformly
    // over a unit sphere, with constraints on cos(theta) and azimuth.
    // These are static so it is only done once, on the first call to
    // the function and is retained for all future calls.
    static double const czmax = 0.5;
    static RandomUnitSphere ranunit(-czmax, czmax);

    // Generate a unit vector in the direction of one track from the origin.
    CLHEP::Hep3Vector vunit(ranunit.shoot());

    // Radius of each layer of the tracker.
    // Access by reference so there is no copying.
    vector<double> const& radius = ctracker->r();

    // Loop over the layers of the CTracker.
    for ( vector<double>::size_type i=0; i<radius.size(); ++i){

      // Path length, s, from origin to intersection with this layer.
      double sinTheta = sqrt( 1. - square(vunit.cosTheta()) );
      double s = radius[i]/sinTheta;

      // Point of intersection between the track and the layer.
      CLHEP::Hep3Vector intersectionPoint = s*vunit;

      // Generate a gaussian distributed pulse height.
      double ph = CLHEP::RandGauss::shoot(10.,2.);
      
      // Keep only hits with large enough pulse height.
      if ( ph > minPulseHeight_ ){
	p->push_back(ToyHit(intersectionPoint,ph));
      }
    }

    // Put the generated hits into the event.
    //event.put(p);

    // At this point p is no longer usable.  If you try to
    // use p it will generate a run time error.

    // There is an alternative syntax for put.
    // Comment out the call to put and uncomment the following lines:
    //
    // Return value is a "read only pointer" to the data product.
    
    edm::OrphanHandle<ToyHitCollection> q = event.put(p);
    edm::LogInfo("Hits") << "Number of hits: " 
			 << q->size();


  }

}


using mu2e::Ex02MakeHits;
DEFINE_FWK_MODULE(Ex02MakeHits);
