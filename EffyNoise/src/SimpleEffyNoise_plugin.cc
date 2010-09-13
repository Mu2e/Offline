// ======================================================================
//
// SimpleEffyNoise_plugin
//
// A class to produce a new StepPointMCCollection, including noise hits
// and excluding hits lost due to inefficiency.
//
// This works on StepPointMC's, and is only intended as a stop-gap
// until hit simulations are more advanced.
//
// ======================================================================

// Framework support:
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
using edm::Event;
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"
#include "FWCore/Services/interface/TFileService.h"

// Mu2e support:
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
using namespace mu2e;

//CLHEP support:
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Vector/ThreeVector.h"
using CLHEP::Hep3Vector;

//Root support:
#include "TH1F.h"

// C++ support:
#include <iostream>
#include <string>
using namespace std;


// ======================================================================


namespace mu2e {

  class SimpleEffyNoise
    : public edm::EDProducer
  {
  public:
    explicit SimpleEffyNoise( edm::ParameterSet const& pset )
    : _diagLevel      ( pset.getUntrackedParameter<int>("diagLevel",0) )
    , _noiseRate      ( pset.getParameter<double>("noiseRate") )
    , _hitIneffy      ( pset.getParameter<double>("hitIneffy") )
    , _messageCategory( "EffyNoise" )
    , _hHitsLostRate  ( 0 )
    , _hHitsLost      ( 0 )
    , _hNoiseHits     ( 0 )
    , _hNoiseHits_Z   ( 0 )
    , _effyFlat       ( createEngine(-1)   )  // use default seed
    , _noiseFlat      ( _effyFlat.engine() )  // share engine
    , _noisePoisson   ( _effyFlat.engine() )  // share engine
    {
      // A place holder.
      produces<StepPointMCCollection>();
    }

    virtual ~SimpleEffyNoise()  { }

    virtual void beginJob( edm::EventSetup const& );

    void produce( edm::Event& e, edm::EventSetup const& );

  private:
    int _diagLevel;     // diagnostics level
    double _noiseRate;  // the per-channel-per-event noise rate
    double _hitIneffy;  // the single hit inefficiency

    // A category for the error logger.
    const std::string _messageCategory;

    TH1F* _hHitsLostRate;// histogram of hits lost
    TH1F* _hHitsLost;   // histogram of hits lost per event
    TH1F* _hNoiseHits;  // histogram of noise hits per event
    TH1F* _hNoiseHits_Z;// histogram of noise hits per event (zoomed)

    CLHEP::RandFlat     _effyFlat;
    CLHEP::RandFlat     _noiseFlat;
    CLHEP::RandPoissonQ _noisePoisson;

  };  // SimpleEffyNoise

  void SimpleEffyNoise::beginJob(edm::EventSetup const& ) {
    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ) {
      edm::Service<edm::TFileService> tfs;

      _hHitsLostRate      = tfs->make<TH1F>( "hHitsLostRate",
          "Hits Lost (0) and Kept (1)",2,-.5,1.5);
      _hHitsLost      = tfs->make<TH1F>( "hHitsLost",
          "Hits Lost per event",50,-.5,49.5);
      _hNoiseHits      = tfs->make<TH1F>( "hNoiseHits",
          "Noise Hits per Event",100,-5,995);
      _hNoiseHits_Z      = tfs->make<TH1F>( "hNoiseHits_Z",
          "Noise Hits per Event",100,-.5,99.5);
    }

  }  // beginJob()

  void SimpleEffyNoise::produce(edm::Event& event, edm::EventSetup const&) {

  std::cout<<"\n\n IN SimpleEffyNoise::produced ================\n\n"<<std::endl;
   // This function has two main sections:
   //  1) apply inefficiency to existing hits
   //  2) add noise hits


    // A container to hold the output hits.
    auto_ptr<StepPointMCCollection> newPoints(new StepPointMCCollection);

    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");
    static const string collectionName("tracker");

    // Ask the event to give us a handle to the requested hits.
    edm::Handle<StepPointMCCollection> points;
    event.getByLabel(creatorName,collectionName,points);

    // Product Id of the input points.
    edm::ProductID const& id( points.id() );

    // first, apply inefficiency
    for ( int ih = 0, sz = points->size(); ih != sz; ++ih ) {
      bool keepPoint = true;
      if (_hitIneffy > 0  &&  _effyFlat.fire() < _hitIneffy)
        keepPoint = false;

      if (keepPoint) {
        const StepPointMC& hit = (*points)[ih];
        newPoints->push_back( StepPointMC(hit) );
        if (_diagLevel>0)  _hHitsLostRate->Fill(1);
      } else {
        if (_diagLevel>0)  _hHitsLostRate->Fill(0);
      }

    } // original hit ih

    if (_diagLevel>0)
        _hHitsLost->Fill(points->size()-newPoints->size());

    // Add noise
std::cout<<"noise_rate="<<_noiseRate<<std::endl;
    if (_noiseRate>0) {

    // need geometry to know how many straws there are to distribute
    // correctly
       // Most of this is copied from rhbob's 1st version of hough tranform

      GeomHandle<LTracker> ltracker;
      int nstraws = ltracker->getAllStraws().size();

      static int noiseMean = static_cast<int>(nstraws*_noiseRate);
std::cout<<"noiseMean="<<noiseMean<<std::endl;
      int numberOfNoiseHits = _noisePoisson.fire(noiseMean);
std::cout<<"numberOfNoiseHits="<<numberOfNoiseHits<<"\n\n"<<std::endl;
      _hNoiseHits->Fill(numberOfNoiseHits);
      _hNoiseHits_Z->Fill(numberOfNoiseHits);

      //add noise hits
      int trackIDnoise = 2;
      double eDepNoise = 0.;
      double timeNoise = 0.;
      CLHEP::Hep3Vector momentumNoise;

      for (int ia=0; ia != numberOfNoiseHits; ++ia) {
        int istraw = static_cast<int>(nstraws*_noiseFlat.fire());
        Straw const& straw = ltracker->getStraw( StrawIndex(istraw) );
        CLHEP::Hep3Vector mid = straw.getMidPoint(); //leftover from Bob's
              // initial HoughTransform work - may need to change.
        CLHEP::Hep3Vector w   = straw.getDirection();

        // Safe dummy value.
        double stepLength(1.0);

        newPoints->push_back(StepPointMC(trackIDnoise,istraw,eDepNoise,timeNoise,0,mid,momentumNoise,stepLength));
      } // noise hit ia

    } // non-zero noise rate

    // save new points and return
    event.put(newPoints);

  } // produce()

}  // namespace mu2e


// ======================================================================


using mu2e::SimpleEffyNoise;
DEFINE_FWK_MODULE(SimpleEffyNoise);

