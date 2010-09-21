//
// An EDProducer Module that checks radiative pi decays
//
// $Id: RPC_plugin.cc,v 1.5 2010/09/21 00:02:36 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/09/21 00:02:36 $
//
// Original author R. Bernstein
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <set>
#include <utility>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "DataFormats/Common/interface/Handle.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"


// Mu2e includes.
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeneralUtilities/inc/pow.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "TrackerGeom/inc/StrawIndex.hh"


//CLHEP includes
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Randomize.h"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::RandPoisson;
using namespace mu2e;
namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 

  class RPC : public edm::EDAnalyzer{
  public:
    RPC(edm::ParameterSet const& pset):

      //
      // Run time parameters
      _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
      _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
      _minimumEnergy(pset.getParameter<double>("minimumEnergy")),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _messageCategory("ToyHitInfo"){}

    virtual ~RPC() {}

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    virtual void beginRun(edm::Run const &r, 
                          edm::EventSetup const& eSetup );

    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lblock, 
                                      edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);


  private:

  

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
 

    // Number of events analyzed.
    int _nAnalyzed;


    // A category for the error logger.
    const std::string _messageCategory;


    TH1D* _piCaptureConvertedElectronMomentum;
    TH1D* _piCaptureConvertedElectronMomentumSignal;
    TH1D* _piCaptureConvertedPositronMomentum;
    TH1D* _piCaptureConvertedPositronMomentumSignal;

    TH1D* _piCaptureConvertedElectronCosTheta;
    TH1D* _piCaptureConvertedElectronCosThetaSignal;
    TH1D* _piCaptureConvertedPositronCosTheta;
    TH1D* _piCaptureConvertedPositronCosThetaSignal;

    TH1D* _conversionAsymmetry;

    TH1D* _piCaptureConvertedElectronHitTracker;
    TH1D* _piCaptureConvertedElectronSignalHitTracker;

    TH1D* _numberOfHitStraws;
    TH1D* _momentumEnteringTracker;
    TH1D* _momentumLostBeforeTracker;

    void RPC::bookEventHistos(double const elow, double const ehigh);
    void RPC::fillEventHistos();
  };


  void RPC::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    //    edm::Service<edm::TFileService> tfs;

  }

  void RPC::endJob(){
  }



  void RPC::beginRun(edm::Run const& run,
		     edm::EventSetup const& eSetup ){
  }

  void RPC::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
				 edm::EventSetup const&){
  }


  void RPC::analyze(const edm::Event& event, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    if (ncalls >= 81000){cout << "ncalls = " << ncalls << endl;} //assert(2==1);


    // Book histogram on the first call regardless
    double elow = 103.5;
    double ehigh = 105.;

    if ( ncalls == 1) bookEventHistos(elow,ehigh);

    edm::LogVerbatim log(_messageCategory);
    log << "RPC event #: " 
        << event.id();

    // 
    // start looking through SimParticles
    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    edm::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    //
    // get handle to hits:  with _trackerStepPoints, if there are hits they belong to the tracker
    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);


    double photonEnergy = 0.;
    double electronEnergy = 0.;
    double positronEnergy = 0.;
    bool electronHitTracker = false;
    bool hitEnoughStraws = false;
    CLHEP::Hep3Vector  momentumAtEntranceToTracker;

    if (haveSimPart){
      for (uint32_t ithPart = 0; ithPart < simParticles->size(); ++ithPart){
        SimParticle const& sim = simParticles->at(ithPart);
	PhysicalVolumeInfo const& startVol = volumes->at(sim.startVolumeIndex());
	//
	// is this the initial photon?
	if (sim.parentId() == -1){
	  if (sim.pdgId() != PDGCode::gamma) {
	    //
	    // this can't happen if we're studying RPCs so throw and die
	    throw cms::Exception("GEOM")
	      << "RPC with a parent not a photon, but parent PDG code is "
	      << sim.pdgId();
	  } else
	    {	CLHEP::HepLorentzVector photonMomentum = sim.startMomentum();
	      photonEnergy = photonMomentum.e();
	    }

	}
	//	cout << " volumename = " << startVol.name() << endl;
	//
	// check three things:  (1) the mother is the original photon, (2) you're an e+ or e-, and (3) the photon converts in the foil
	if ( sim.parentId() == 0 && startVol.name() == "TargetFoil_" ){
	  //	if ( sim.parentId() == 0 && startVol.name() == "ToyDSCoil" ){
	  if (sim.pdgId() == PDGCode::e_minus) {

            //
            // before proceeding I want to be sure the electron made it to the tracker.  Did it hit?  I also need to check the hit came from
            // the generated electron
            std::set<StrawIndex> hitStraws;

            for( size_t i=0; i<hits->size(); ++i ){

              const StepPointMC& hit = (*hits)[i];

              //step point mc associated with generated electrons     
              int trackId = hit.trackId();
              //
              // now I have the track Id of this hit.  
              SimParticle const& simParticleHit = simParticles->at(trackId);
              //
              // is this simParticle associated with the hit the same as the one I started with?
              if (sim.id() == simParticleHit.id()){
                electronHitTracker = true;
                //
                // let's plot how many different straws this electron hit. Set doesn't let me insert duplicates
                hitStraws.insert(hit.strawIndex());

                //
                // the momentum at the first hit of the tracker is the momentum entering the tracker
                if (i==0){
                  momentumAtEntranceToTracker = hit.momentum();
                  //                  cout << "momentum at entrance to tracker = " << momentumAtEntranceToTracker.mag() << endl;
                  _momentumEnteringTracker->Fill(momentumAtEntranceToTracker.mag());
                }

                /*
                  cout << "strawPair: " 
                  << event.id() << " "
                  << i << " "
                  << hit.strawIndex() << " "
                  << strawPair.second << " "
                  << hitStraws.size() << " "
                  << numberOfHitStraws
                  << endl;
                */
              }

            } // end loop over hits.
            if (electronHitTracker) {
              _numberOfHitStraws->Fill(static_cast<double>(hitStraws.size()));
          

              const CLHEP::HepLorentzVector& electronMomentum = sim.startMomentum();
              electronEnergy = electronMomentum.e();
              double momentum = sqrt(pow(electronMomentum.e(),2) - pow(electronMomentum.invariantMass(),2));
              _piCaptureConvertedElectronMomentum->Fill(momentum);
              //
              // crude attempt to model acceptance.  This 0.15 is about a 40% cut.  But then not clear how to use 0.8 cut from MECO-140R.  Still,
              // an interesting number.
              if (abs(electronMomentum.cosTheta()) <= 0.15 && momentum > elow && momentum < ehigh)
                {_piCaptureConvertedElectronMomentumSignal->Fill(momentum);}

              //
              // plot momentum for electrons that hit the tracker
                _piCaptureConvertedElectronHitTracker->Fill(momentum);
                if (momentum > elow && momentum < ehigh)_piCaptureConvertedElectronSignalHitTracker->Fill(momentum);
              _piCaptureConvertedElectronCosTheta->Fill( electronMomentum.cosTheta() );
              if (momentum > elow && momentum < ehigh) _piCaptureConvertedElectronCosThetaSignal->Fill(electronMomentum.cosTheta() );
              //
              //what's the energy loss between birth and entrance?
              _momentumLostBeforeTracker->Fill(momentum - momentumAtEntranceToTracker.mag());


	  
            }
            if (sim.pdgId() == PDGCode::e_plus) {
              CLHEP::HepLorentzVector electronMomentum = sim.startMomentum();
              positronEnergy = electronMomentum.e();
              double momentum = sqrt(pow(electronMomentum.e(),2) - pow(electronMomentum.invariantMass(),2));
              _piCaptureConvertedPositronMomentum->Fill(momentum);
              if(momentum > elow && momentum < ehigh) _piCaptureConvertedPositronMomentumSignal->Fill(momentum);
              _piCaptureConvertedPositronCosTheta->Fill(electronMomentum.cosTheta() );
              if (momentum > elow && momentum < ehigh) _piCaptureConvertedPositronCosThetaSignal->Fill(electronMomentum.cosTheta() );
            }
          }
        }
        // 
        // check geant4's photon conversion.  if there was no conversion you won't see anything so check
        if (simParticles->size() >= 3 && electronEnergy > 0. && positronEnergy > 0. && photonEnergy > 0.){
          _conversionAsymmetry->Fill( abs(electronEnergy - positronEnergy)/ (electronEnergy + positronEnergy) );
        }
      }
    }
  }
  void RPC::bookEventHistos(double const elow, double const ehigh)
  {        
    //    cout << "booking histos" << endl; assert(2==1);
    edm::Service<edm::TFileService> tfs;
    _piCaptureConvertedElectronMomentum = 
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentum",
                       "Pi Capture Converted Electron Momentum", 200, 0., 200.);
    _piCaptureConvertedElectronMomentumSignal = 
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumSignal",
                       "Pi Capture Converted Electron MomentumSignal", 15, elow, ehigh);

    _piCaptureConvertedElectronHitTracker = 
      tfs->make<TH1D>( "piCaptureConvertedElectronHitTracker",
                       "Pi Capture Converted Electron Momentum and Hit Tracker", 200, 0., 200.);


    _piCaptureConvertedElectronSignalHitTracker = 
      tfs->make<TH1D>( "piCaptureConvertedElectronSignalHitTracker",
                       "Pi Capture Converted Electron Momentum and Hit Tracker", 15, elow, ehigh);



    _piCaptureConvertedPositronMomentum = 
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentum",
                       "Pi Capture Converted Positron Momentum", 200, 0., 200.);
    _piCaptureConvertedPositronMomentumSignal = 
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentumSignal",
                       "Pi Capture Converted Positron MomentumSignal", 20, elow, ehigh);

    _piCaptureConvertedElectronCosTheta = 
      tfs->make<TH1D>( "piCaptureConvertedElectronCosTheta",
                       "Pi Capture Converted Electron CosTheta", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaSignal = 
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaSignal",
                       "Pi Capture Converted Electron CosThetaSignal", 200, -1., 1.);
    _piCaptureConvertedPositronCosTheta = 
      tfs->make<TH1D>( "piCaptureConvertedPositronCosTheta",
                       "Pi Capture Converted Positron CosTheta", 200, -1., 1.);
    _piCaptureConvertedPositronCosThetaSignal = 
      tfs->make<TH1D>( "piCaptureConvertedPositronCosThetaSignal",
                       "Pi Capture Converted Positron CosThetaSignal", 200, -1., 1.);

    _conversionAsymmetry = 
      tfs->make<TH1D>("conversionAsymmetry","electron-positron/ephoton", 100,0.,1.);

    _numberOfHitStraws =
      tfs->make<TH1D>("numberOfHitStraws"," Number Of Hit Straws", 100,0.,100.);

    _momentumEnteringTracker =
      tfs->make<TH1D>("momentumEnteringTracker","Electron Momentum Entering Tracker", 40,100.,110.);
    _momentumLostBeforeTracker =
      tfs->make<TH1D>("momentumLostBeforeTracker","Electron Momentum Lost Before Tracker", 25,0.,5.);


  } //bookEventHistos


} // namespace mu2e


DEFINE_FWK_MODULE(RPC);
