/*----------------------------------------------------------------------

  A producer module that makes a collection of overly simplified "hits"
  and adds them to the event.

  $Id: Ex02MakeHits_plugin.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
  $Author: kutschke $
  $Date: 2009/09/30 22:57:47 $
   
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
#include "TargetGeom/inc/Target.hh"
#include "CTrackerGeom/inc/CTracker.hh"

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
    
    virtual void beginRun(edm::Run &r, edm::EventSetup const& eSetup );
    
    static void fillDescription(edm::ParameterSetDescription& iDesc,
                                string const& moduleLabel) {
      iDesc.setAllowAnything();
    }
    
    // Generate one track from the origin.
    void genTrack();

    // A uniform random number generator.
    // Should use the one provided by the framework - not working today.
    // Returns a uniform random number on [xl,xh] or (xl,xh)?
    double myRandom( double xl=0, double xh=1.);
    
    // An approximate gaussian random number generator.
    // Should use the one provided by the framework - not working today.
    double myRandomGauss(double mean, double sigma);

  private:

    // Parameters from run time configuration.
    double minPulseHeight_;
    
    // Properties of the generated track.
    // Filled in genTrack(). Used in produce().
    double phi_;   // Azimuth
    double cphi_;  // cos(Azimuth)
    double sphi_;  // sin(Azimuth)
    double ct_;    // cos(polar angle)
    double st_;    // sin(polar angle)
    double cot_;   // cotan(polar angle)
    
  };

  // At begin run time, update the geometry.
  // In the future this will change to a geometryChanged() method.
  void Ex02MakeHits::beginRun( edm::Run &r, edm::EventSetup const& eSetup ){
    
    // Get access to the geometry system.
    edm::Service<GeometryService> geom;
    
    // Check that the detector parts that we want are present.
  }
  
  void
  Ex02MakeHits::produce(edm::Event& e, edm::EventSetup const&) {
    GeomHandle<Target> target;
    GeomHandle<CTracker> ctracker;

    // Extract the event number from the event.
    int eventNumber = e.id().event();
    
    // Generate the direction of one track from the origin.
    genTrack();

    // Make the collection to hold the hits that we will make.
    auto_ptr<ToyHitCollection> p(new ToyHitCollection);

    // Radius of each layer of the tracker.
    // Access by const reference so there is no copying.
    const vector<double>& radius = ctracker->r();

    // Loop over the layers of the CTracker.
    for ( int i=0; i<radius.size(); ++i){

      // Compute the point at which the straight track intercepts this layer.
      double r  = radius[i];
      double x  = r*cphi_;
      double y  = r*sphi_;
      double z  = r*cot_;

      // Generate a gaussian distributed pulse height.
      double ph = myRandomGauss(10.,2.);
      
      // Keep only hits with large enough pulse height.
      if ( ph > minPulseHeight_ ){
	p->push_back(ToyHit(x,y,z,ph));
      }
    }

    // How many hits were generated?
    int nhits(p->size());

    // Put the hits into the event.  You may no longer use p.
    // Return value is a "read only pointer" to the hits.
    edm::OrphanHandle<ToyHitCollection> q = e.put(p);

    // A silly example to show that you can use q.
    if ( q->size() != nhits ){
      edm::LogError("InconsistentData")
	<< "Number of hits changed after I put them into the event.  "
	<< "  Before: " << nhits
	<< "  After: " << q->size();
    }

    // If you do not need q, there is an simpler syntax to put p into the event.
    //e.put(p);

  }

  // Remove this when we have a real random number generator.
  // Flat random value on the interval (xl,xh].
  double Ex02MakeHits::myRandom(double xl, double xh){
    double x = drand48();
    double l = (xh-xl);
    if ( l <= 0. ){
      throw edm::Exception(edm::errors::UnimplementedFeature,
			   "Hacked random number generator is broken.");      
    }
    return xl + x*l;
  }

  // Remove this when we have a real random number generator.
  // Approximate method to generate a random value drawn from a 
  // gaussian with specified mean and sigma
  double Ex02MakeHits::myRandomGauss(double mean, double sigma){
    double sum(0.);
    for ( int i=0; i<12; ++i){
      sum += myRandom(-0.5,0.5);
    }
    return mean+sum*sigma;
  }

  // Generate the direction of one track, flat in phase space but
  // cut far from the z axis.
  void Ex02MakeHits::genTrack(){
    phi_  = myRandom( 0., 2.*M_PI);
    ct_   = myRandom( -0.5, 0.5);
    cphi_ = cos(phi_);
    sphi_ = sin(phi_);
    st_   = sqrt( 1. - ct_*ct_);

    // I was careful to define ct in a range that is safe for this.
    cot_  = ct_/st_;

  }

}


using mu2e::Ex02MakeHits;
DEFINE_FWK_MODULE(Ex02MakeHits);
